//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XolotlUserObject.h"
#include <dlfcn.h>
//MOOSE includes

#ifdef __APPLE__
#define ISSTANDALONE true
#elif __linux__
#define ISSTANDALONE false
#endif

registerMooseObject("MooseApp", XolotlUserObject);

template <>
InputParameters
validParams<XolotlUserObject>()
{
  InputParameters params = validParams<NodalUserObject>();
  params.addClassDescription("Executes the external application solving for a simple diffusion equation using Finite Difference method. The mesh parameters should be copied from GeneratedMesh and pasted here.");
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addCoupledVar("variable_gb", "The name of the variable that passes the grain boundary data");
  params.addParam<std::string>("library_path_name",
                               "default",
                               "Name with the path for the dynamic library to load");
  params.addParam<std::string>("XolotlInput_path_name",
                               "default",
                               "Name with the path for the Xolotl input file");

  return params;
}

XolotlUserObject::XolotlUserObject(const InputParameters & parameters)
  : NodalUserObject(parameters),
  MooseVariableInterface<Real>(this,
                               false,
                               "variable",
                               Moose::VarKindType::VAR_ANY,
                               Moose::VarFieldType::VAR_FIELD_STANDARD),
  _ext_coord(_mesh.nNodes()),
  _ext_data(_mesh.nNodes()),
  _var(*mooseVariable()),
  _u(_var.dofValues()),
  _v(coupledValue("variable_gb")),
  _ext_lib_path_name(getParam<std::string>("library_path_name")),
  _xolotl_input_path_name(getParam<std::string>("XolotlInput_path_name"))
{
  std::cout<<"Initialization of XolotlUserObject"<<std::endl;

  /* Finite Difference mesh parameters for the External App.
     using MOOSE mesh parameters assigned in GeneratedMesh.
     Only works for uniform QUAD4 mesh.
     Assuming that MOOSE mesh and the mesh of the External App. has the same geometry.
     initExtCoords() generates a finite difference mesh that contains an ordered set of contiguous memory address
       --> Depends on how the External App. solves the finite difference problem.
       --> In this example, the External App. calculates the spatial second derivative as (f[i-1] + 2.*f[i] + f[i+1])/dx/dx;
       --> Thus the contiguous memory is required.
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &_moose_rank);

  // Loading a fake Xolotl command line arguments
  _argv[0] = new char[_parameterFile.length()+1];
  strcpy(_argv[0], _parameterFile.c_str());
  _argv[1] = new char[_xolotl_input_path_name.length() + 1];
  strcpy(_argv[1], _xolotl_input_path_name.c_str());
  _argv[2] = 0; // null-terminate the array

  // _ext_coord = initExtCoords();
  // print_ext_coord(_ext_coord);

  // Preparing for the external library handle
  _ext_lib_handle = dlopen(_ext_lib_path_name.c_str(), RTLD_LAZY);
  if (!_ext_lib_handle) {
    std::cerr << "Cannot open library: " << dlerror() << '\n';
  }
  dlerror();

  // Loading Xolotl solver pointer
  typedef XolotlDLinterface* create_t();
  typedef void destroy_t(XolotlDLinterface*);

  create_t* create_interface = (create_t*) dlsym(_ext_lib_handle, "create");
  const char* dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol create: " << dlsym_error << '\n';
  }

  destroy_t* destroy_interface = (destroy_t*) dlsym(_ext_lib_handle, "destroy");
  dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol destroy: " << dlsym_error << '\n';
  }

  // create an instance of the class
  // XolotlDLinterface* interface = create_interface();
  _xolotl_interface = create_interface();

  // Loading Xolotl options from the Xolotl input file
  // _xolotl_options.readParams(_argv);

  _xolotl_solver = _xolotl_interface->initializeXolotl(_argc, _argv, MPI_COMM_WORLD, ISSTANDALONE);

  // Get Xolotl grid information
  _xolotl_interface->getXolotlGlobalGridInfo(&_xolotl_dim, &_xolotl_nx, &_xolotl_ny, &_xolotl_nz,
  &_xolotl_dx, &_xolotl_dy, &_xolotl_dz, &_xolotl_lx, &_xolotl_ly, &_xolotl_lz, _argc, _argv);

  _xolotl_xc = build_xolotl_axis(_xolotl_nx, _xolotl_dx);
  _xolotl_yc = build_xolotl_axis(_xolotl_ny, _xolotl_dy);
  _xolotl_zc = build_xolotl_axis(_xolotl_nz, _xolotl_dz);
  for (int i = 0; i < _xolotl_ny; i++) {
    printf("%d\t%lf\n", i, _xolotl_yc[i]);
  }


  // Print out the loaded Xolotl grid parameters
  print_mesh_params();

  _xolotl_local_index_bounds = init_xolotl_local_index_table(6);
  fillout_xolotl_local_index_table(_xolotl_local_index_bounds, 6);
  // print_xolotl_local_index_table(_xolotl_local_index_bounds, 6);

  //  Get local array index bounds
  int xs, xm, Mx, ys, ym, My, zs, zm, Mz;
  _xolotl_interface->getLocalCoordinates(_xolotl_solver, &xs, &xm, &Mx, &ys, &ym, &My, &zs, &zm, &Mz);
  _xolotl_xi_lb = xs;
  _xolotl_xi_ub = xs + xm;
  _xolotl_yi_lb = ys;
  _xolotl_yi_ub = ys + ym;
  _xolotl_zi_lb = zs;
  _xolotl_zi_ub = zs + zm;

  // Syncing the time stepping
  Real dtime = 1.1e-20;
  if (_dt > 1e-20) dtime = _dt;
  _xolotl_interface->setTimes(_xolotl_solver, _t, dtime);
  _xolotl_LocalXeRate = _xolotl_interface->getLocalXeRate(_xolotl_solver);  // Buffer initialization; data values do not matter here
  _xolotl_LocalConc = _xolotl_interface->getLocalXeConc(_xolotl_solver);  // Buffer initialization; data values do not matter here


  // Printing out the nodeID through the domain
  // std::vector<dof_id_type> nodelist;
  // nodelist = _mesh.getNodeList();
  // std::cout<<"MPIrank = "<<_moose_rank<<", "<<"Size of node list: "<<nodelist.size()<<std::endl;
  // for (auto i = 0; i < nodelist.size(); i++)
  // {
  //   Node * nd = _mesh.queryNodePtr(i); //This call retrives a Node pointer corresponding to the Node ID provided.
  //   unsigned int nodeID = nd->id();
  //   std::cout<<"Node loop; id: "<<nodeID<<std::endl;
  // }

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // std::cout << "constructor: MPIrank = " << myrank <<std::endl;
  // printf("constructor: MPIrank = %d\n", myrank);
}

void
XolotlUserObject::initialize()
{
  // The parameters of the External Application will be initialized here
  // This member function is called once at a MOOSE time step (maybe able to be specified by using execute_on parameter in the input file)

  // Syncing the time stepping
  _xolotl_interface->setTimes(_xolotl_solver, _t, _dt);

  int xolotl_local_nx = _xolotl_LocalXeRate->size();
  int xolotl_local_ny = _xolotl_LocalXeRate->at(0).size();
  int xolotl_local_nz = _xolotl_LocalXeRate->at(0)[0].size();

  // Console output of the data and its coordinate
  for (int k = 0; k < xolotl_local_nz; k++){
    for (int j = 0; j < xolotl_local_ny; j++){
      for (int i = 0; i < xolotl_local_nx; i++){
        double localRate = _xolotl_LocalXeRate->at(i)[j][k];
        if (fabs(localRate) > 1e-20) {
          std::cout << i<<", "<<j<<", "<<k<<", "<<localRate << std::endl;
        }
      }
    }
  }

  int xs, xm, Mx, ys, ym, My, zs, zm, Mz;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  _xolotl_interface->getLocalCoordinates(_xolotl_solver, &xs, &xm, &Mx, &ys, &ym, &My, &zs, &zm, &Mz);
  std::cout<<"rank, (xs, xm), (ys, ym), (zs, zm)"<<std::endl;
  printf("%d, (%d, %d), (%d, %d), (%d, %d)\n",myrank, xs, xm, ys, ym, zs, zm);

  // std::cout<<

  // std::cout << "initialize: MPIrank = " << myrank <<std::endl;
  // printf("initialize: MPIrank = %d\n", myrank);
}

void
XolotlUserObject::execute()
{
  // Data transfer from MOOSE to External App.
  // This member function is called at every MOOSE node points.
  // _ext_data is initialized by copying the values from the assigned AuxVariable

  // unsigned int ii = map_MOOSE2Xolotl(*_current_node);

  // Initializing _ext_data with the AuxVariable value; _u[0]
  // _ext_data[ii] = _u[0];
  // Get MPI rank
  // int myrank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //
  // unsigned int nodeID = _current_node->id();
  //std::cout << "execute: MPIrank, nodeID = " << myrank <<", "<<nodeID << std::endl;
  // printf("execute: MPIrank = %d, nodeID = %d\n", myrank, nodeID);

}

void
XolotlUserObject::finalize()
{
  _xolotl_interface->solveXolotl(_xolotl_solver);
  _xolotl_LocalXeRate = _xolotl_interface->getLocalXeRate(_xolotl_solver); // Bringing the Xe rate data (within GB)
  _xolotl_LocalConc = _xolotl_interface->getLocalXeConc(_xolotl_solver);  // Bringing the Xe concentration data (within bulk)

  // Updating _ext_data by executing the External App.
  // This member function is called once at a time step.
  // The external app will be instanciated and update the _ext_data here

  // Normal run
  // Temporary buffers for the External App.
  // double *mdp_data;
  // mdp_data = new double[mdp.numdata];
  // double *mdp_x;
  // mdp_x = new double[mdp.numdata];
  //
  // //Copying the coordinates and data to the temporary buffers for the External App.
  // for (unsigned int i = mdp.iLower; i <= mdp.iUpper; i++) {
  //   mdp_x[i] = _ext_coord[i-mdp.nghost](0); // Copying x-coordinates; beware that the index of mdp_x has been shifted as mdp.nghost
  //   mdp_data[i] = _ext_data[i-mdp.nghost]; // Copying the data; beware that the index of mdp_data has also been shifted as mdp.nghost
  // }

  //Printing the data just for checking; before the update
  // mdp_output_data_console(mdp_data, mdp.iLower, mdp.iUpper, mdp_x, mdp.time_start);

  // Updating mdp_data
  // mdp_Diffu_module(mdp_data, mdp_x, mdp);

  // //Printing the data just for checking; after the update
  // mdp_output_data_console(mdp_data, mdp.iLower, mdp.iUpper, mdp_x, mdp.time_end);
  //
  // //Updating _ext_data with mdp_data updated
  // for (unsigned int i = mdp.iLower; i <= mdp.iUpper; i++) {
  //   _ext_data[i-mdp.nghost] = mdp_data[i]; // Copying the data; beware that the index of mdp_data has also been shifted as mdp.nghost
  // }
  //
  // // Deallocating the temporary buffers
  // delete[] mdp_data;
  // delete[] mdp_x;

  // _xolotl_interface->finalizeXolotl(_xolotl_solver, false);
  // _xolotl_solver.reset();
  // Deallocating the external library handle
  // dlclose(_ext_lib_handle);
  // int myrank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // std::cout << "finalize: MPIrank = " << myrank <<std::endl;
  // printf("finalize: MPIrank = %d\n", myrank);
}

void
XolotlUserObject::threadJoin(const UserObject & /*y*/)
{
  //I don't know what to do with this yet.
}

Real
XolotlUserObject::calc_spatial_value() const
{
  // Data transfer from External App. to MOOSE
  // This member function is called as per spatialValue() called, which is called at every node points.
  // The corresponding ext_data to the current node will be retunred by this function

  // unsigned int ii = map_MOOSE2Ext(*_current_node);

  // Returning _ext_data values to pass to the AuxKernel that will modify the AuxVariable values
  // return _ext_data[ii];
  // int myrank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //std::cout << "calc_spatial_value: MPIrank = " << myrank <<s

  int rank, i, j, k;
  map_MOOSE2Xolotl(&rank, &i, &j, &k, *_current_node);
  double localRate = _xolotl_LocalXeRate->at(i)[j][k];
  double concentration = _xolotl_LocalConc->at(i)[j][k];

  // return localRate;
  return concentration;
}

// Real
// XolotlUserObject::getMinInDimension(unsigned int component) const
// {
//   switch (component)
//   {
//     case 0:
//       return _xmin;
//     case 1:
//       return _dim > 1 ? _ymin : 0;
//     case 2:
//       return _dim > 2 ? _zmin : 0;
//     default:
//       mooseError("Invalid component");
//   }
// }
//
// Real
// XolotlUserObject::getMaxInDimension(unsigned int component) const
// {
//   switch (component)
//   {
//     case 0:
//       return _xmax;
//     case 1:
//       return _dim > 1 ? _ymax : 0;
//     case 2:
//       return _dim > 2 ? _zmax : 0;
//     default:
//       mooseError("Invalid component");
//   }
// }
//
// std::vector<libMesh::Point>
// XolotlUserObject::initExtCoords() const
// {
//   std::vector<libMesh::Point> xyz(_mesh.nNodes());
//   print_mesh_params();
//   // Contiguous memory chunk for the Ext. app. (i-outermost order)
//   for (unsigned int i = 0; i < nNode_x; i++) {
//     for (unsigned int j = 0; j < nNode_y; j++) {
//       for (unsigned int k = 0; k < nNode_z; k++) {
//         unsigned int ii = nNode_y * nNode_z * i + nNode_z * j + k;
//         libMesh::Point xyz_build={0.0};
//         xyz_build(0) = _xmin + i * _xolotl_dx;
//         xyz_build(1) = _ymin + j * _xolotl_dy;
//         xyz_build(2) = _zmin + k * _xolotl_dz;
//         xyz[ii] = xyz_build;
//         std::cout<<"x,y,z = "<<xyz[ii](0)<<","<<xyz[ii](1)<<","<<xyz[ii](2)<<","<<std::endl;
//       }
//     }
//   }
//   return xyz;
// }
//
void
XolotlUserObject::map_MOOSE2Xolotl(int *rankreturn, int *ireturn, int *jreturn, int *kreturn, const Node & MOOSEnode) const
{
  double xn = MOOSEnode(0);
  double yn = MOOSEnode(1);
  double zn = MOOSEnode(2);
  int iglob = 0, jglob = 0, kglob = 0; //Global grid indices of Xolotl, which is the nearest grid point from the given MOOSE node.
  double distx, disty, distz;
  int isave, jsave, ksave;
  int maxsize = max3int(_xolotl_nx, _xolotl_ny, _xolotl_nz);
  //Serching for iglob
  distx = _xolotl_lx; //initialize with the maximum distance
  disty = _xolotl_ly; //initialize with the maximum distance
  distz = _xolotl_lz; //initialize with the maximum distance
  isave = -1;
  jsave = -1;
  ksave = -1;
  for (int i = 0; i < maxsize; i++) {
    //Searching for iglob
    if (i < _xolotl_nx) {
      double d = fabs(_xolotl_xc[i] - xn);
      if (d <= distx) {
        distx = d;
        isave = i;
      }
    }
    if (i < _xolotl_ny) {
      double d = fabs(_xolotl_yc[i] - yn);
      if (d <= disty) {
        disty = d;
        jsave = i;
      }
    }
    if (i < _xolotl_nz) {
      double d = fabs(_xolotl_zc[i] - zn);
      if (d <= distz) {
        distz = d;
        ksave = i;
      }
    }
  }

  iglob = isave;
  jglob = jsave;
  kglob = ksave;

  // Determing Xolotl MPI rank of the processor where iglob, jglob, and kglob are belong
  int nprocess;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
  int ranksave = -1;
  for (int i = 0; i < nprocess; i++) {
    int ixL = _xolotl_local_index_bounds[i][0];
    int ixU = _xolotl_local_index_bounds[i][1];
    int iyL = _xolotl_local_index_bounds[i][2];
    int iyU = _xolotl_local_index_bounds[i][3];
    int izL = _xolotl_local_index_bounds[i][4];
    int izU = _xolotl_local_index_bounds[i][5];
    if ( iglob >= ixL && iglob <= ixU
      && jglob >= iyL && jglob <= iyU
      && kglob >= izL && kglob <= izU )
      {
        ranksave = i;
        break;
      }
  }

  // Converting iglob, jglob, and kglob into the local indices
  int iloc = iglob - _xolotl_local_index_bounds[ranksave][0];
  int jloc = jglob - _xolotl_local_index_bounds[ranksave][2];
  int kloc = kglob - _xolotl_local_index_bounds[ranksave][4];

  *rankreturn = ranksave;
  if (ranksave == _moose_rank) {
    *ireturn = iloc;
    *jreturn = jloc;
    *kreturn = kloc;
  }else{
    std::cout<<"Parallel transfer is not supported yet."<<std::endl;
    return;
  }
}

void
XolotlUserObject::print_mesh_params() const
{
  std::cout<<"_xolotl_dim = "<<_xolotl_dim<<std::endl;
  std::cout<<"_xolotl_lx = "<<_xolotl_lx<<std::endl;
  std::cout<<"_xolotl_ly = "<<_xolotl_ly<<std::endl;
  std::cout<<"_xolotl_lz = "<<_xolotl_lz<<std::endl;
  std::cout<<"_xolotl_nx, _xolotl_ny, _xolotl_nz = "<<_xolotl_nx<<","<<_xolotl_ny<<","<<_xolotl_nz<<std::endl;
  std::cout<<"_xolotl_dx, _xolotl_dy, _xolotl_dz = "<<_xolotl_dx<<","<<_xolotl_dy<<","<<_xolotl_dz<<std::endl;
}
//
// void
// XolotlUserObject::print_ext_coord(std::vector<libMesh::Point> a) const
// {
//   for(unsigned int i = 0; i < a.size() ; i++)
//   {
//     std::cout<<"x,y,z = "<<a[i](0)<<","<<a[i](1)<<","<<a[i](2)<<std::endl;
//   }
// }

int**
XolotlUserObject::init_xolotl_local_index_table(int ncolumn) const
{
  int nprocess;
  int **table;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
  table = new int*[nprocess];
  for (int i = 0; i < nprocess; i++) {
    table[i] = new int[ncolumn];
  }
  for (int i = 0; i < nprocess; i++) {
    for (int j = 0; j < ncolumn; j++) {
      table[i][j] = 0;
    }
  }
  return table;
}

void
XolotlUserObject::print_xolotl_local_index_table(int **table, int ncol) const
{
  int nrow;
  MPI_Comm_size(MPI_COMM_WORLD, &nrow);
  for (int nrank = 0; nrank < nrow; nrank++) {
    printf("Mrank = %d, ( Xrank = %d, ", _moose_rank, nrank);
    for (int j = 0; j < ncol; j++) {
      printf("%d%s",table[nrank][j], j == ncol - 1 ? ")\n":", ");
    }
  }
}

void
XolotlUserObject::fillout_xolotl_local_index_table(int **table, int ncol) const
{
  int nrow;
  MPI_Comm_size(MPI_COMM_WORLD, &nrow);
  int xs, xm, Mx, ys, ym, My, zs, zm, Mz;
  _xolotl_interface->getLocalCoordinates(_xolotl_solver, &xs, &xm, &Mx, &ys, &ym, &My, &zs, &zm, &Mz);
  table[_moose_rank][0] = xs;
  table[_moose_rank][1] = xs + xm;
  table[_moose_rank][2] = ys;
  table[_moose_rank][3] = ys + ym;
  table[_moose_rank][4] = zs;
  table[_moose_rank][5] = zs + zm;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int globalsum = 0;
      MPI_Allreduce(&table[i][j], &globalsum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
      table[i][j] = globalsum;
    }
  }

}

double*
XolotlUserObject::build_xolotl_axis(int nsize, double dl) const
{
  double *tmpbuff;
  tmpbuff = new double[nsize];
  for (int i = 0; i < nsize; i++) {
    tmpbuff[i] = (double) i * dl;
  }
  return tmpbuff;
}

int
XolotlUserObject::max3int(int a, int b, int c) const
{
  int result = a;
  if (b >= result) {
    result = b;
  }
  if (c >= result) {
    result = c;
  }
  return result;
}
