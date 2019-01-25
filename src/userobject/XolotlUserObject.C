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
  // _ext_coord(_mesh.nNodes()),
  // _ext_data(_mesh.nNodes()),
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
  _xolotl_interface = create_interface();

  _xolotl_solver = _xolotl_interface->initializeXolotl(_argc, _argv, MPI_COMM_WORLD, ISSTANDALONE);

  // Get Xolotl grid information
  _xolotl_interface->getXolotlGlobalGridInfo(&_xolotl_dim, &_xolotl_nx, &_xolotl_ny, &_xolotl_nz,
  &_xolotl_dx, &_xolotl_dy, &_xolotl_dz, &_xolotl_lx, &_xolotl_ly, &_xolotl_lz, _argc, _argv);

  _xolotl_xc = build_xolotl_axis(_xolotl_nx, _xolotl_dx);
  _xolotl_yc = build_xolotl_axis(_xolotl_ny, _xolotl_dy);
  _xolotl_zc = build_xolotl_axis(_xolotl_nz, _xolotl_dz);

  // Print out the loaded Xolotl grid parameters
  print_mesh_params();

  _xolotl_local_index_table = init_xolotl_local_index_table(6);
  fillout_xolotl_local_index_table(_xolotl_local_index_table, 6);
  print_xolotl_local_index_table(_xolotl_local_index_table, 6);

  _xolotl_xi_lb = _xolotl_local_index_table[_moose_rank][0];
  _xolotl_xi_ub = _xolotl_xi_lb + _xolotl_local_index_table[_moose_rank][1] - 1;
  _xolotl_yi_lb = _xolotl_local_index_table[_moose_rank][2];
  _xolotl_yi_ub = _xolotl_yi_lb + _xolotl_local_index_table[_moose_rank][3] - 1;
  _xolotl_zi_lb = _xolotl_local_index_table[_moose_rank][4];
  _xolotl_zi_ub = _xolotl_zi_lb + _xolotl_local_index_table[_moose_rank][5] - 1;
  _xolotl_localNx = _xolotl_local_index_table[_moose_rank][1];
  _xolotl_localNy = _xolotl_local_index_table[_moose_rank][3];
  _xolotl_localNz = _xolotl_local_index_table[_moose_rank][5];

  // Syncing the time stepping
  Real dtime = 1.1e-20;
  if (_dt > 1e-20) dtime = _dt;
  _xolotl_interface->setTimes(_xolotl_solver, _t, dtime);
  _xolotl_XeRate = vectorized_xolotl_XeRate(_xolotl_solver);
  _xolotl_XeConc = vectorized_xolotl_XeConc(_xolotl_solver);
  _xolotl_GlobalXeRate = allocate_xolotlGlobalData();
  _xolotl_GlobalXeConc = allocate_xolotlGlobalData();

  // Make the local buffer visible for all processors in MPI_COMM_WORLD
  // MPI_Win_create(&_xolotl_XeConc[0], _xolotl_localNx * _xolotl_localNy * _xolotl_localNz * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &_win);

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

  // Simple MPI_get example

}

void
XolotlUserObject::initialize()
{
  // The parameters of the External Application will be initialized here
  // This member function is called once at a MOOSE time step (maybe able to be specified by using execute_on parameter in the input file)

  // Syncing the time stepping
  _xolotl_interface->setTimes(_xolotl_solver, _t, _dt);

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
  _xolotl_XeRate = vectorized_xolotl_XeRate(_xolotl_solver);
  _xolotl_XeConc = vectorized_xolotl_XeConc(_xolotl_solver);
  localFill_xolotlGlobalXeRate(_xolotl_GlobalXeRate, _xolotl_solver);
  globalFill_xolotlGlobalData(_xolotl_GlobalXeRate);
  localFill_xolotlGlobalXeConc(_xolotl_GlobalXeConc, _xolotl_solver);
  globalFill_xolotlGlobalData(_xolotl_GlobalXeConc);

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

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  Real result;
  int recvRank, i, j, k;
  map_MOOSE2Xolotl(&recvRank, &i, &j, &k, *_current_node);
  if (myrank == recvRank) { // MPI communication will not happen.
    result = _xolotl_XeConc[ii(i,j,k)];
  }else{  //MPI communication will happen.
    double tmp = -99999.9;
    // MPI_Win win;
    // MPI_Win_create(_xolotl_XeConc, _xolotl_localNx * _xolotl_localNy * _xolotl_localNz * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    // // if (myrank == 2) printf("myrank=%d, recvRank=%d, i,j,k = %d, %d, %d\n",myrank,recvRank,i,j,k);
    // // MPI_Win_fence(0, win); // One-sided communication begin
    // MPI_Win_lock(MPI_LOCK_SHARED, recvRank, 0, win);
    // int err = MPI_Get(&tmp, 1, MPI_DOUBLE, recvRank, iiR(recvRank,i,j,k), 1, MPI_DOUBLE, win);
    // // // MPI_Win_fence(0, win); // One-sided communication end
    // MPI_Win_unlock(recvRank, win);
    // MPI_Win_free(&win);
    // if (err != MPI_SUCCESS) {
    //   if (err == MPI_ERR_ARG) {
    //     printf("MPI_Get:: invalid argument, rank = %d\n)", myrank);
    //   }
    //   if (err == MPI_ERR_COUNT) {
    //     printf("MPI_Get:: invalid count, rank = %d\n)", myrank);
    //   }
    //   if (err == MPI_ERR_RANK) {
    //     printf("MPI_Get:: invalid rank, rank = %d\n)", myrank);
    //   }
    //   if (err == MPI_ERR_TYPE) {
    //     printf("MPI_Get:: invalid type, rank = %d\n)", myrank);
    //   }
    //   if (err == MPI_ERR_WIN) {
    //     printf("MPI_Get:: invalid window, rank = %d\n)", myrank);
    //   } else {
    //     printf("MPI_Get:: unknown error, rank = %d\n", myrank);
    //   }
    // }else{
    //   // if (myrank == 2)
    //     printf("MPI_Get:: Success!, myrank=%d, recvRank=%d, i,j,k = %d, %d, %d\n", myrank,recvRank,i,j,k);
    // }
    result = tmp;
  }


  // double tmp = -9999.9;
  // // MPI_Win_fence(MPI_MODE_NOPRECEDE, _win);
  // MPI_Win_fence(0, win);
  // int err = MPI_Get(&tmp, 1, MPI_DOUBLE, recvRank, iiR(recvRank,i,j,k), 1, MPI_DOUBLE, win);
  // MPI_Win_fence(0, win);
  // if (err != MPI_SUCCESS) {
  //   if (err == MPI_ERR_ARG) {
  //     printf("MPI_Get:: invalid argument, rank = %d\n)", myrank);
  //   }
  //   if (err == MPI_ERR_COUNT) {
  //     printf("MPI_Get:: invalid count, rank = %d\n)", myrank);
  //   }
  //   if (err == MPI_ERR_RANK) {
  //     printf("MPI_Get:: invalid rank, rank = %d\n)", myrank);
  //   }
  //   if (err == MPI_ERR_TYPE) {
  //     printf("MPI_Get:: invalid type, rank = %d\n)", myrank);
  //   }
  //   if (err == MPI_ERR_WIN) {
  //     printf("MPI_Get:: invalid window, rank = %d\n)", myrank);
  //   } else {
  //     printf("MPI_Get:: unknown error, rank = %d\n", myrank);
  //   }
  // }else{
  //   printf("MPI_Get:: Success!, tmp = %lf, myrank=%d, recvRank=%d, i,j,k = %d, %d, %d, iiR = %d\n",tmp, myrank,recvRank,i,j,k, iiR(recvRank, i, j, k));
  // }
  // result = tmp;
  // // MPI_Win_fence(MPI_MODE_NOSUCCEED, _win);
  // MPI_Win_free(&win);
  // return result;
  return (double)myrank;
}

Real
XolotlUserObject::calc_spatial_value_glob() const
{
  int i, j, k;
  map_MOOSE2XolotlGlob(&i, &j, &k, *_current_node);
  return _xolotl_GlobalXeConc[iiGlob(i,j,k)];
  // return _moose_rank;
}

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
    int ixL = _xolotl_local_index_table[i][0];
    int ixU = ixL + _xolotl_local_index_table[i][1] - 1;
    int iyL = _xolotl_local_index_table[i][2];
    int iyU = iyL + _xolotl_local_index_table[i][3] - 1;
    int izL = _xolotl_local_index_table[i][4];
    int izU = izL + _xolotl_local_index_table[i][5] - 1;
    if ( iglob >= ixL && iglob <= ixU
      && jglob >= iyL && jglob <= iyU
      && kglob >= izL && kglob <= izU )
      {
        ranksave = i;
        break;
      }
  }

  // Converting iglob, jglob, and kglob into the local indices
  int iloc = iglob - _xolotl_local_index_table[ranksave][0];
  int jloc = jglob - _xolotl_local_index_table[ranksave][2];
  int kloc = kglob - _xolotl_local_index_table[ranksave][4];

  *rankreturn = ranksave;
  *ireturn = iloc;
  *jreturn = jloc;
  *kreturn = kloc;
}

void
XolotlUserObject::map_MOOSE2XolotlGlob(int *ireturn, int *jreturn, int *kreturn, const Node & MOOSEnode) const
{
  double xn = MOOSEnode(0);
  double yn = MOOSEnode(1);
  double zn = MOOSEnode(2);
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

  *ireturn = isave;
  *jreturn = jsave;
  *kreturn = ksave;
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
  // table[_moose_rank][1] = xm;
  table[_moose_rank][1] = xm < 1 ? 1 : xm;
  table[_moose_rank][2] = ys;
  // table[_moose_rank][3] = ym;
  table[_moose_rank][3] = ym < 1 ? 1 : ym;
  table[_moose_rank][4] = zs;
  // table[_moose_rank][5] = zm;
  table[_moose_rank][5] = zm < 1 ? 1 : zm;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int globalsum = 0;
      MPI_Allreduce(&table[i][j], &globalsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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

int**
XolotlUserObject::init_xolotl_rankpair() const
{
  int nprocess;
  int ncolumn = 2;
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
XolotlUserObject::fillout_xolotl_rankpair(int **table, int myrank, int recvRank) const
{
  int nrow;
  MPI_Comm_size(MPI_COMM_WORLD, &nrow);
  int ncol = 2;
  printf("fillout_xolotl_rankpair:: _moose_rank=%d, recvRank=%d\n",myrank, recvRank);
  table[recvRank][0] = myrank; //Filling out the sendRank
  table[myrank][1] = recvRank; //Filling out the recvRank
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int globalsum = 0;
      MPI_Allreduce(&table[i][j], &globalsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      table[i][j] = globalsum;
    }
  }

}

void
XolotlUserObject::print_xolotl_rankpair(int **table) const
{
  int nrow;
  MPI_Comm_size(MPI_COMM_WORLD, &nrow);
  for (int nrank = 0; nrank < nrow; nrank++) {
    printf("%d, %d, %d\n", table[nrank][0], nrank, table[nrank][1]);
  }
}

double*
XolotlUserObject::vectorized_xolotl_XeRate(std::shared_ptr<xolotlSolver::PetscSolver> solver) const
{
  std::vector<std::vector<std::vector<double> > > * buff = _xolotl_interface->getLocalXeRate(solver); // Bringing the Xe rate data (within GB)

  int nx = _xolotl_localNx;
  int ny = _xolotl_localNy;
  int nz = _xolotl_localNz;

  double *vecReturn;
  vecReturn = new double[nx*ny*nz];
  for (int k = 0; k < nz; k++){
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        vecReturn[ii(i,j,k)] = buff->at(i)[j][k];
      }
    }
  }

  return vecReturn;
}

double*
XolotlUserObject::vectorized_xolotl_XeConc(std::shared_ptr<xolotlSolver::PetscSolver> solver) const
{
  std::vector<std::vector<std::vector<double> > > * buff = _xolotl_interface->getLocalXeConc(solver); // Bringing the Xe rate data (within GB)
  int nx = _xolotl_localNx;
  int ny = _xolotl_localNy;
  int nz = _xolotl_localNz;
  double *vecReturn;
  vecReturn = new double[nx*ny*nz];
  for (int k = 0; k < nz; k++){
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        vecReturn[ii(i,j,k)] = buff->at(i)[j][k];
      }
    }
  }

  return vecReturn;
}

double*
XolotlUserObject::allocate_xolotlGlobalData() const
{
  int nsize = _xolotl_nx * _xolotl_ny * _xolotl_nz;
  double *arr = new double[nsize];
  for (int i = 0; i < nsize; i++){
    arr[i] = 0.0;
  }

  return arr;
}

void
XolotlUserObject::localFill_xolotlGlobalXeRate(double *arr, std::shared_ptr<xolotlSolver::PetscSolver> solver) const
{
  std::vector<std::vector<std::vector<double> > > * buff = _xolotl_interface->getLocalXeRate(solver); // Bringing the Xe rate data (within GB)
  // for (int k = _xolotl_zi_lb; k <= _xolotl_zi_ub; k++) {
  //   for (int j = _xolotl_yi_lb; j <= _xolotl_yi_ub; j++) {
  //     for (int i = _xolotl_xi_lb; i <= _xolotl_xi_ub; i++) {
  //       int iloc = i - _xolotl_xi_lb;
  //       int jloc = j - _xolotl_yi_lb;
  //       int kloc = k - _xolotl_zi_lb;
  //       arr[iiGlob(i,j,k)] = buff->at(iloc)[jloc][kloc];
  //     }
  //   }
  // }
  int nsize = _xolotl_nx * _xolotl_ny * _xolotl_nz;
  for (int ii = 0; ii < nsize; ii++){
    arr[ii] = 0.0;
  }
  for (int k = 0; k < _xolotl_localNz; k++){
    for (int j = 0; j < _xolotl_localNy; j++){
      for (int i = 0; i < _xolotl_localNx; i++){
        int iGlob = i +_xolotl_xi_lb;
        int jGlob = j +_xolotl_yi_lb;
        int kGlob = k +_xolotl_zi_lb;
        arr[iiGlob(iGlob,jGlob,kGlob)] = buff->at(i)[j][k];
      }
    }
  }
}

void
XolotlUserObject::localFill_xolotlGlobalXeConc(double *arr, std::shared_ptr<xolotlSolver::PetscSolver> solver) const
{
  std::vector<std::vector<std::vector<double> > > * buff = _xolotl_interface->getLocalXeConc(solver); // Bringing the Xe rate data (within GB)
  // for (int k = _xolotl_zi_lb; k <= _xolotl_zi_ub; k++) {
  //   for (int j = _xolotl_yi_lb; j <= _xolotl_yi_ub; j++) {
  //     for (int i = _xolotl_xi_lb; i <= _xolotl_xi_ub; i++) {
  //       int iloc = i - _xolotl_xi_lb;
  //       int jloc = j - _xolotl_yi_lb;
  //       int kloc = k - _xolotl_zi_lb;
  //       arr[iiGlob(i,j,k)] = buff->at(iloc)[jloc][kloc];
  //     }
  //   }
  // }
  // Cleaning up the buffer
  int nsize = _xolotl_nx * _xolotl_ny * _xolotl_nz;
  for (int ii = 0; ii < nsize; ii++){
    arr[ii] = 0.0;
  }
  for (int k = 0; k < _xolotl_localNz; k++){
    for (int j = 0; j < _xolotl_localNy; j++){
      for (int i = 0; i < _xolotl_localNx; i++){
        int iGlob = i +_xolotl_xi_lb;
        int jGlob = j +_xolotl_yi_lb;
        int kGlob = k +_xolotl_zi_lb;
        // if (_moose_rank == 1) {
        //   arr[iiGlob(iGlob,jGlob,kGlob)] = buff->at(i)[j][k];
        // }
        arr[iiGlob(iGlob,jGlob,kGlob)] = buff->at(i)[j][k];
        // arr[iiGlob(iGlob,jGlob,kGlob)] = 1.0;
        // arr[iiGlob(iGlob,jGlob,kGlob)] = _moose_rank;
      }
    }
  }
}

void
XolotlUserObject::globalFill_xolotlGlobalData(double *arr) const
{
  int nsize = _xolotl_nx * _xolotl_ny * _xolotl_nz;
  for (int i = 0; i < nsize; i++) {
    double globalsum = 0.0;
    MPI_Allreduce(&arr[i], &globalsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    arr[i] = globalsum;
  }
  // MPI_Barrier(MPI_COMM_WORLD);
}

int
XolotlUserObject::iiGlob(int i, int j, int k) const
{
  return _xolotl_ny * _xolotl_nz * i + _xolotl_nz * j + k;
  // return i + _xolotl_nx * j + _xolotl_nx * _xolotl_ny * k;
}

int
XolotlUserObject::ii(int i, int j, int k) const
{
  return _xolotl_localNy * _xolotl_localNz * i + _xolotl_localNz * j + k;
  // return i + _xolotl_localNx * j + _xolotl_localNx * _xolotl_localNy * k;
}

int
XolotlUserObject::iiR(int rank, int i, int j, int k) const
{
  int nx = _xolotl_local_index_table[rank][1];
  int ny = _xolotl_local_index_table[rank][3];
  int nz = _xolotl_local_index_table[rank][5];
  return ny * nz * i + nz * j + k;
  // return i + nx * j + nx * ny * k;
}
