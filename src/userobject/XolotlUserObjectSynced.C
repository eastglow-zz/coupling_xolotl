//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XolotlUserObjectSynced.h"
#include <dlfcn.h>
//MOOSE includes

#ifdef __APPLE__
#define ISSTANDALONE true
#elif __linux__
#define ISSTANDALONE false
#endif

registerMooseObject("MooseApp", XolotlUserObjectSynced);

template <>
InputParameters
validParams<XolotlUserObjectSynced>()
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
  params.addParam<Real>("gb_marker_threshold",
                               0.95,
                               "Threshold to identify the non-matrix region (GB or bubble) when (bnds(x,y,z) <= gb_marker_threshold), where bnds = 1 within the matrix region, 0 within the bubble region, and 0.5 <= bnds <= 1.0 within grain boundaries.");

  return params;
}

XolotlUserObjectSynced::XolotlUserObjectSynced(const InputParameters & parameters)
  : NodalUserObject(parameters),
  MooseVariableInterface<Real>(this,
                               false,
                               "variable",
                               Moose::VarKindType::VAR_ANY,
                               Moose::VarFieldType::VAR_FIELD_STANDARD),
  _var(*mooseVariable()),
  _u(_var.dofValues()),
  _v(coupledValue("variable_gb")),
  _ext_lib_path_name(getParam<std::string>("library_path_name")),
  _xolotl_input_path_name(getParam<std::string>("XolotlInput_path_name")),
  _matrix_marker_thres(getParam<Real>("gb_marker_threshold"))
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
  _xolotl_interface->getXolotlGlobalGridInfo(&_xolotl_dim, &_xolotl_regulargrid,
  &_xolotl_nx, &_xolotl_ny, &_xolotl_nz,
  &_xolotl_dx, &_xolotl_dy, &_xolotl_dz, &_xolotl_lx, &_xolotl_ly, &_xolotl_lz, _argc, _argv);

  _xolotl_xc = build_xolotl_axis(_xolotl_nx, _xolotl_dx);
  _xolotl_yc = build_xolotl_axis(_xolotl_ny, _xolotl_dy);
  _xolotl_zc = build_xolotl_axis(_xolotl_nz, _xolotl_dz);
  if (!_xolotl_regulargrid) {
    _xolotl_xcNR = _xolotl_interface->getXolotlXgrid(_xolotl_solver);
    _xolotl_lxNR = _xolotl_xcNR[_xolotl_nx-1];
  }

  // Print out the loaded Xolotl grid parameters
  // print_mesh_params();

  _xolotl_local_index_table = init_xolotl_local_index_table(6);
  fillout_xolotl_local_index_table(_xolotl_local_index_table, 6);
  // print_xolotl_local_index_table(_xolotl_local_index_table, 6);

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

  _xolotl_XeRate = allocate_xolotlLocalData();
  _xolotl_GlobalXeRate = allocate_xolotlGlobalData();
}

void
XolotlUserObjectSynced::initialize()
{
  // The parameters of the External Application will be initialized here
  // This member function is called once at a MOOSE time step (maybe able to be specified by using execute_on parameter in the input file)

  // Syncing the time stepping
  // printf("rank = %d, _dt_old = %lf, _t = %lf, _dt = %lf\n", _moose_rank, _dt_old, _t, _dt);
  _xolotl_interface->setTimes(_xolotl_solver, _t, _dt);
  _xolotl_interface->initGBLocation(_xolotl_solver);
}

void
XolotlUserObjectSynced::execute()
{
  if (_v[0] <= _matrix_marker_thres) {
    int i, j, k;
    map_MOOSE2XolotlGlob(&i, &j, &k, *_current_node);
    _GBListLocal.push_back(std::make_tuple(i,j,k));
  }

}

void
XolotlUserObjectSynced::finalize()
{
  //Get XeRateOld
  double *XeRateOld = vectorized_xolotl_XeRate(_xolotl_solver);

  _GBList = get_GlobalGBList(_GBListLocal);
  _xolotl_interface->setGBLocations(_xolotl_solver, _GBList);

  _xolotl_interface->solveXolotl(_xolotl_solver);

  //Get XeRateNew
  double *XeRateNew = vectorized_xolotl_XeRate(_xolotl_solver);

  double *Rate = computeTimeDerivativeLocal(XeRateNew, XeRateOld, _dt);

  localFill_xolotlGlobalData(_xolotl_GlobalXeRate, Rate);

  std::vector<std::tuple<int, int, int>> GBListLocalClean;
  _GBListLocal = GBListLocalClean;
}

void
XolotlUserObjectSynced::threadJoin(const UserObject & /*y*/)
{
  //I don't know what to do with this yet.
}

Real
XolotlUserObjectSynced::calc_spatial_value() const
{
  // Data transfer from External App. to MOOSE
  // This member function is called as per spatialValue() called, which is called at every node points.
  // The corresponding ext_data to the current node will be retunred by this function

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // Real result;
  // int recvRank, i, j, k;
  // map_MOOSE2Xolotl(&recvRank, &i, &j, &k, *_current_node);
  // if (myrank == recvRank) { // MPI communication will not happen.
  //   result = _xolotl_XeConc[ii(i,j,k)];
  // }else{  //MPI communication will happen.
  //   double tmp = -99999.9;
  //   result = tmp;
  // }

  return (double)myrank;
}

Real
XolotlUserObjectSynced::calc_spatial_value_glob() const
{
  int i, j, k;
  map_MOOSE2XolotlGlob(&i, &j, &k, *_current_node);
  return _xolotl_GlobalXeRate[iiGlob(i,j,k)];
}

void
XolotlUserObjectSynced::map_MOOSE2Xolotl(int *rankreturn, int *ireturn, int *jreturn, int *kreturn, const Node & MOOSEnode) const
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
XolotlUserObjectSynced::map_MOOSE2XolotlGlob(int *ireturn, int *jreturn, int *kreturn, const Node & MOOSEnode) const
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
      double x = _xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i];
      double d = fabs(x - xn);
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
XolotlUserObjectSynced::print_mesh_params() const
{
  std::cout<<"_xolotl_dim = "<<_xolotl_dim<<std::endl;
  if (_xolotl_regulargrid) {
    std::cout<<"_xolotl_lx = "<<_xolotl_lx<<std::endl;
  }else{
    std::cout<<"_xolotl_lxNR = "<<_xolotl_lxNR<<std::endl;
  }
  std::cout<<"_xolotl_ly = "<<_xolotl_ly<<std::endl;
  std::cout<<"_xolotl_lz = "<<_xolotl_lz<<std::endl;
  std::cout<<"_xolotl_nx, _xolotl_ny, _xolotl_nz = "<<_xolotl_nx<<","<<_xolotl_ny<<","<<_xolotl_nz<<std::endl;
  if (_xolotl_regulargrid) {
    std::cout<<"_xolotl_dx, _xolotl_dy, _xolotl_dz = "<<_xolotl_dx<<","<<_xolotl_dy<<","<<_xolotl_dz<<std::endl;
  }else{
    std::cout<<"_xolotl_dy, _xolotl_dz = "<<_xolotl_dy<<","<<_xolotl_dz<<std::endl;
  }
}

int**
XolotlUserObjectSynced::init_xolotl_local_index_table(int ncolumn) const
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
XolotlUserObjectSynced::print_xolotl_local_index_table(int **table, int ncol) const
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
XolotlUserObjectSynced::fillout_xolotl_local_index_table(int **table, int ncol) const
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

// double*
std::vector<double>
XolotlUserObjectSynced::build_xolotl_axis(int nsize, double dl) const
{
  // double *tmpbuff;
  // tmpbuff = new double[nsize];
  std::vector<double> tmpbuff;
  for (int i = 0; i < nsize; i++) {
    // tmpbuff[i] = (double) i * dl;
    tmpbuff.push_back((double) i * dl);
  }
  return tmpbuff;
}

int
XolotlUserObjectSynced::max3int(int a, int b, int c) const
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


double*
XolotlUserObjectSynced::allocate_xolotlLocalData() const
{
  int nsize = _xolotl_localNx * _xolotl_localNy * _xolotl_localNz;
  double *arr = new double[nsize];
  for (int i = 0; i < nsize; i++){
    arr[i] = 0.0;
  }

  return arr;
}

double*
XolotlUserObjectSynced::vectorized_xolotl_XeRate(std::shared_ptr<xolotlSolver::PetscSolver> solver) const
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
XolotlUserObjectSynced::vectorized_xolotl_XeConc(std::shared_ptr<xolotlSolver::PetscSolver> solver) const
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
XolotlUserObjectSynced::computeTimeDerivativeLocal(double *buff_new, double *buff_old, double timeInterval) const
{
  int localsize = _xolotl_localNx * _xolotl_localNy * _xolotl_localNz;
  double *toReturn = new double[localsize];
  for (int i = 0; i < localsize; i++) {
    if (timeInterval > 0) {
      toReturn[i] = (buff_new[i] - buff_old[i])/timeInterval;
    }else{
      toReturn[i] = 0.0;
      printf("XolotlUserObject::computeTimeDerivativeLocal(): divide by zero or negative!\n");
    }
  }
  return toReturn;
}

void
XolotlUserObjectSynced::superposeArrayLocal(double *ans, double *arr1, double *arr2) const
{
  int localsize = _xolotl_localNx * _xolotl_localNy * _xolotl_localNz;
  for (int i = 0; i < localsize; i++) {
    ans[i] = arr1[i] + arr2[i];
  }
}

double*
XolotlUserObjectSynced::allocate_xolotlGlobalData() const
{
  int nsize = _xolotl_nx * _xolotl_ny * _xolotl_nz;
  double *arr = new double[nsize];
  for (int i = 0; i < nsize; i++){
    arr[i] = 0.0;
  }

  return arr;
}

void
XolotlUserObjectSynced::localFill_xolotlGlobalData(double *arr, double *arrLocal) const
{
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
        arr[iiGlob(iGlob,jGlob,kGlob)] = arrLocal[ii(i,j,k)];
      }
    }
  }
}

void
XolotlUserObjectSynced::localFill_xolotlGlobalXeRate(double *arr, std::shared_ptr<xolotlSolver::PetscSolver> solver) const
{
  std::vector<std::vector<std::vector<double> > > * buff = _xolotl_interface->getLocalXeRate(solver); // Bringing the Xe rate data (within GB)

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
XolotlUserObjectSynced::localFill_xolotlGlobalXeConc(double *arr, std::shared_ptr<xolotlSolver::PetscSolver> solver) const
{
  std::vector<std::vector<std::vector<double> > > * buff = _xolotl_interface->getLocalXeConc(solver); // Bringing the Xe rate data (within GB)

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
        arr[iiGlob(iGlob,jGlob,kGlob)] = buff->at(i)[j][k];
      }
    }
  }
}

void
XolotlUserObjectSynced::globalFill_xolotlGlobalData(double *arr) const
{
  int nsize = _xolotl_nx * _xolotl_ny * _xolotl_nz;
  for (int i = 0; i < nsize; i++) {
    double globalsum = 0.0;
    MPI_Allreduce(&arr[i], &globalsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    arr[i] = globalsum;
  }
}

std::vector<std::tuple<int, int, int>>
XolotlUserObjectSynced::get_GlobalGBList(std::vector<std::tuple<int, int, int>> gbLocal) const
{
  int localsize = gbLocal.size();
  int globalsize = 0;

  //Make the size table
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int *sizes = new int[nprocs];
  for (int i = 0; i < nprocs; i ++) {
    sizes[i] = 0;
  }
  sizes[_moose_rank] = localsize;
  for (int i = 0; i < nprocs; i++ ) {
    int tmp = 0;
    MPI_Allreduce(&sizes[i], &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    sizes[i] = tmp;
  }

  //Get the total GB coordinate list size
  MPI_Allreduce(&localsize, &globalsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  //Initializing the global GB coordinate list
  int *igb_global = new int[globalsize];
  int *jgb_global = new int[globalsize];
  int *kgb_global = new int[globalsize];
  for (int i = 0; i < globalsize; i++) {
    igb_global[i] = 0.0;
    jgb_global[i] = 0.0;
    kgb_global[i] = 0.0;
  }

  // Calc. displacement to save the GB coordinate list
  int idisp = 0;
  for (int i = 0; i < _moose_rank && _moose_rank > 0; i++) {
    idisp += sizes[i];
  }

  if (localsize > 0) {
    int *igb = new int[localsize];
    int *jgb = new int[localsize];
    int *kgb = new int[localsize];

    // Get local GB coordinates
    for (int i = 0; i < localsize; i++) {
      igb[i] = std::get<0>(gbLocal[i]);
      jgb[i] = std::get<1>(gbLocal[i]);
      kgb[i] = std::get<2>(gbLocal[i]);
    }

    //Fill the local GB coord. data in the global array
    for (int i = idisp; i < idisp + localsize; i++) {
      igb_global[i] = igb[i-idisp];
      jgb_global[i] = jgb[i-idisp];
      kgb_global[i] = kgb[i-idisp];
    }
  }

  //Merge the local GB coords into the global GB coord array
  for (int i = 0; i < globalsize; i++) {
    int tmpi = 0, tmpj = 0, tmpk = 0;
    MPI_Allreduce(&igb_global[i], &tmpi, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    igb_global[i] = tmpi;
    MPI_Allreduce(&jgb_global[i], &tmpj, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    jgb_global[i] = tmpj;
    MPI_Allreduce(&kgb_global[i], &tmpk, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    kgb_global[i] = tmpk;
  }

  //Make coordinates into tuple to return
  std::vector<std::tuple<int, int, int>> gbListGlobal;
  for (int i = 0; i < globalsize; i++) {
    gbListGlobal.push_back(std::make_tuple(igb_global[i], jgb_global[i], kgb_global[i]));
  }

  return gbListGlobal;
}

int
XolotlUserObjectSynced::iiGlob(int i, int j, int k) const
{
  // return _xolotl_ny * _xolotl_nz * i + _xolotl_nz * j + k;
  return i + _xolotl_nx * j + _xolotl_nx * _xolotl_ny * k;
}

int
XolotlUserObjectSynced::ii(int i, int j, int k) const
{
  // return _xolotl_localNy * _xolotl_localNz * i + _xolotl_localNz * j + k;
  return i + _xolotl_localNx * j + _xolotl_localNx * _xolotl_localNy * k;
}

int
XolotlUserObjectSynced::iiR(int rank, int i, int j, int k) const
{
  int nx = _xolotl_local_index_table[rank][1];
  int ny = _xolotl_local_index_table[rank][3];
  int nz = _xolotl_local_index_table[rank][5];
  // return ny * nz * i + nz * j + k;
  return i + nx * j + nx * ny * k;
}
