//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XolotlMesh.h"

#include "libmesh/mesh_generation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/unstructured_mesh.h"

// C++ includes
#include <cmath> // provides round, not std::round (see http://www.cplusplus.com/reference/cmath/round/)
#include <dlfcn.h>

#ifdef __APPLE__
#define ISSTANDALONE true
#elif __linux__
#define ISSTANDALONE false
#endif

registerMooseObject("MooseApp", XolotlMesh);

template <>
InputParameters
validParams<XolotlMesh>()
{
  InputParameters params = validParams<MooseMesh>();

  MooseEnum elem_types(
      "EDGE EDGE2 EDGE3 EDGE4 QUAD QUAD4 QUAD8 QUAD9 TRI3 TRI6 HEX HEX8 HEX20 HEX27 TET4 TET10 "
      "PRISM6 PRISM15 PRISM18 PYRAMID5 PYRAMID13 PYRAMID14"); // no default

  MooseEnum dims("1=1 2 3");
  params.addRequiredParam<MooseEnum>(
      "dim", dims, "The dimension of the mesh to be generated"); // Make this parameter required

  // params.addParam<unsigned int>("nx", 1, "Number of elements in the X direction");
  // params.addParam<unsigned int>("ny", 1, "Number of elements in the Y direction");
  // params.addParam<unsigned int>("nz", 1, "Number of elements in the Z direction");
  // params.addParam<Real>("xmin", 0.0, "Lower X Coordinate of the generated mesh");
  // params.addParam<Real>("ymin", 0.0, "Lower Y Coordinate of the generated mesh");
  // params.addParam<Real>("zmin", 0.0, "Lower Z Coordinate of the generated mesh");
  // params.addParam<Real>("xmax", 1.0, "Upper X Coordinate of the generated mesh");
  // params.addParam<Real>("ymax", 1.0, "Upper Y Coordinate of the generated mesh");
  // params.addParam<Real>("zmax", 1.0, "Upper Z Coordinate of the generated mesh");
  params.addParam<MooseEnum>("elem_type",
                             elem_types,
                             "The type of element from libMesh to "
                             "generate (default: linear element for "
                             "requested dimension)");
  params.addParam<bool>(
      "gauss_lobatto_grid",
      false,
      "Grade mesh into boundaries according to Gauss-Lobatto quadrature spacing.");
  params.addRangeCheckedParam<Real>(
      "bias_x",
      1.,
      "bias_x>=0.5 & bias_x<=2",
      "The amount by which to grow (or shrink) the cells in the x-direction.");
  params.addRangeCheckedParam<Real>(
      "bias_y",
      1.,
      "bias_y>=0.5 & bias_y<=2",
      "The amount by which to grow (or shrink) the cells in the y-direction.");
  params.addRangeCheckedParam<Real>(
      "bias_z",
      1.,
      "bias_z>=0.5 & bias_z<=2",
      "The amount by which to grow (or shrink) the cells in the z-direction.");
  params.addParam<std::string>("library_path_name",
      "default",
      "Name with the path for the dynamic library to load");
  params.addParam<std::string>("XolotlInput_path_name",
      "default",
      "Name with the path for the Xolotl input file");

  // params.addParamNamesToGroup("dim", "Main");

  params.addClassDescription(
      "Create a line, square, or cube mesh with uniformly spaced or biased elements.");
  return params;
}

XolotlMesh::XolotlMesh(const InputParameters & parameters)
  : MooseMesh(parameters),
    _dim(getParam<MooseEnum>("dim")),
    // _nx(getParam<unsigned int>("nx")),
    // _ny(getParam<unsigned int>("ny")),
    // _nz(getParam<unsigned int>("nz")),
    // _xmin(getParam<Real>("xmin")),
    // _xmax(getParam<Real>("xmax")),
    // _ymin(getParam<Real>("ymin")),
    // _ymax(getParam<Real>("ymax")),
    // _zmin(getParam<Real>("zmin")),
    // _zmax(getParam<Real>("zmax")),
    _gauss_lobatto_grid(getParam<bool>("gauss_lobatto_grid")),
    _bias_x(getParam<Real>("bias_x")),
    _bias_y(getParam<Real>("bias_y")),
    _bias_z(getParam<Real>("bias_z")),
    _ext_lib_path_name(getParam<std::string>("library_path_name")),
    _xolotl_input_path_name(getParam<std::string>("XolotlInput_path_name"))
{
  if (_gauss_lobatto_grid && (_bias_x != 1.0 || _bias_y != 1.0 || _bias_z != 1.0))
    mooseError("Cannot apply both Gauss-Lobatto mesh grading and biasing at the same time.");

  // All generated meshes are regular orthogonal meshes
  _regular_orthogonal_mesh = true;

  /*
    ======= Xolotl related routines below =======
  */
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

  // Get Xolotl grid information
  _xolotl_interface->getXolotlGlobalGridInfo(&_xolotl_dim, &_xolotl_regulargrid,
    &_xolotl_nx, &_xolotl_ny, &_xolotl_nz,
    &_xolotl_dx, &_xolotl_dy, &_xolotl_dz, &_xolotl_lx, &_xolotl_ly, &_xolotl_lz, _argc, _argv);

  _xolotl_xc = build_xolotl_axis(_xolotl_nx, _xolotl_dx);
  _xolotl_yc = build_xolotl_axis(_xolotl_ny, _xolotl_dy);
  _xolotl_zc = build_xolotl_axis(_xolotl_nz, _xolotl_dz);

  _xolotl_solver = _xolotl_interface->initializeXolotl(_argc, _argv, MPI_COMM_WORLD, ISSTANDALONE);
  if (!_xolotl_regulargrid) {
    _xolotl_xcNR = _xolotl_interface->getXolotlXgrid(_xolotl_solver);
  }

  _dim = _xolotl_dim;
  _nx = _xolotl_nx-1; // Note that _nx is the # of elements, and _xolotl_nx is the # of nodes.
  _ny = _xolotl_ny-1;
  _nz = _xolotl_nz-1;
  _xmin = 0.0;
  _xmax = _xolotl_lx;
  _ymin = 0.0;
  _ymax = _xolotl_ly;
  _zmin = 0,0;
  _zmax = _xolotl_lz;

  /*
    ======= Xolotl related routines above =======
  */

}

Real
XolotlMesh::getMinInDimension(unsigned int component) const
{
  switch (component)
  {
    case 0:
      return _xmin;
    case 1:
      return _dim > 1 ? _ymin : 0;
    case 2:
      return _dim > 2 ? _zmin : 0;
    default:
      mooseError("Invalid component");
  }
}

Real
XolotlMesh::getMaxInDimension(unsigned int component) const
{
  switch (component)
  {
    case 0:
      return _xmax;
    case 1:
      return _dim > 1 ? _ymax : 0;
    case 2:
      return _dim > 2 ? _zmax : 0;
    default:
      mooseError("Invalid component");
  }
}

std::unique_ptr<MooseMesh>
XolotlMesh::safeClone() const
{
  return libmesh_make_unique<XolotlMesh>(*this);
}

void
XolotlMesh::buildMesh()
{
  MooseEnum elem_type_enum = getParam<MooseEnum>("elem_type");

  if (!isParamValid("elem_type"))
  {
    // Switching on MooseEnum
    switch (_dim)
    {
      case 1:
        elem_type_enum = "EDGE2";
        break;
      case 2:
        elem_type_enum = "QUAD4";
        break;
      case 3:
        elem_type_enum = "HEX8";
        break;
    }
  }

  ElemType elem_type = Utility::string_to_enum<ElemType>(elem_type_enum);

  // Switching on MooseEnum
  switch (_dim)
  {
    // The build_XYZ mesh generation functions take an
    // UnstructuredMesh& as the first argument, hence the dynamic_cast.
    case 1:
      MeshTools::Generation::build_line(dynamic_cast<UnstructuredMesh &>(getMesh()),
                                        _nx,
                                        _xmin,
                                        _xmax,
                                        elem_type,
                                        _gauss_lobatto_grid);
      break;
    case 2:
      MeshTools::Generation::build_square(dynamic_cast<UnstructuredMesh &>(getMesh()),
                                          _nx,
                                          _ny,
                                          _xmin,
                                          _xmax,
                                          _ymin,
                                          _ymax,
                                          elem_type,
                                          _gauss_lobatto_grid);
      break;
    case 3:
      MeshTools::Generation::build_cube(dynamic_cast<UnstructuredMesh &>(getMesh()),
                                        _nx,
                                        _ny,
                                        _nz,
                                        _xmin,
                                        _xmax,
                                        _ymin,
                                        _ymax,
                                        _zmin,
                                        _zmax,
                                        elem_type,
                                        _gauss_lobatto_grid);
      break;
  }

  // When Xolotl uses non-regular grid
  if (!_xolotl_regulargrid)
  {
    // Reference to the libmesh mesh
    MeshBase & mesh = getMesh();

    // "width" of the mesh in each direction
    Real width[3] = {_xmax - _xmin, _ymax - _ymin, _zmax - _zmin};

    // Min mesh extent in each direction.
    Real mins[3] = {_xmin, _ymin, _zmin};

    // Number of elements in each direction.
    unsigned int nelem[3] = {_nx, _ny, _nz};

    // Loop over the nodes and move them to the desired location
    for (auto & node_ptr : mesh.node_ptr_range())
    {
      Node & node = *node_ptr;

      int i, j, k;
      //Get the global coordinate index of Xolotl grid
      map_MOOSE2XolotlGlob(&i, &j, &k, node);
      Real newcoord[3] = {_xolotl_xcNR[i], _xolotl_yc[j], _xolotl_zc[k]};

      for (unsigned int dir = 0; dir < LIBMESH_DIM; ++dir)
      {
        //Move node to non-regular grid
        node(dir) = mins[dir] + newcoord[dir];
      }
    }
  }
}

std::vector<double>
XolotlMesh::build_xolotl_axis(int nsize, double dl) const
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

void
XolotlMesh::map_MOOSE2XolotlGlob(int *ireturn, int *jreturn, int *kreturn, const Node & MOOSEnode) const
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

int
XolotlMesh::max3int(int a, int b, int c) const
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
