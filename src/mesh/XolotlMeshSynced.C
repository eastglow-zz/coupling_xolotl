//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XolotlMeshSynced.h"

#include "libmesh/mesh_generation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/partitioner.h"
#include "libmesh/metis_csr_graph.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/remote_elem.h"
#include "libmesh/face_quad4.h"
#include "libmesh/cell_hex8.h"

// C++ includes
#include <cmath> // provides round, not std::round (see http://www.cplusplus.com/reference/cmath/round/)
#include <dlfcn.h>

#ifdef __APPLE__
#define ISSTANDALONE true
#elif __linux__
#define ISSTANDALONE false
#endif

registerMooseObject("MooseApp", XolotlMeshSynced);

template <>
InputParameters
validParams<XolotlMeshSynced>()
{
  InputParameters params = validParams<MooseMesh>();

  MooseEnum elem_types("EDGE2  QUAD4  HEX8"); // no default

  MooseEnum dims("1=1 2 3", "2");
  params.addRequiredParam<MooseEnum>(
      "dim", dims, "The dimension of the mesh to be generated"); // Make this parameter required

  params.addParam<MooseEnum>("elem_type",
                             elem_types,
                             "The type of element from libMesh to "
                             "generate (default: linear element for "
                             "requested dimension)");

  params.addParamNamesToGroup("dim", "Main");

  params.addClassDescription(
      "Create a line, square, or cube mesh with using PETSc DMDA.");

  // This mesh is always distributed
  params.set<MooseEnum>("parallel_type") = "DISTRIBUTED";

  // Parameter for the Xolotl file name
  params.addParam<std::string>("library_path_name",
      "default",
      "Name with the path for the dynamic library to load");
  params.addParam<std::string>("XolotlInput_path_name",
      "default",
      "Name with the path for the Xolotl input file");

  return params;
}

XolotlMeshSynced::XolotlMeshSynced(const InputParameters & parameters)
  : MooseMesh(parameters),
    _dim(getParam<MooseEnum>("dim")),
    _ext_lib_path_name(getParam<std::string>("library_path_name")),
    _xolotl_input_path_name(getParam<std::string>("XolotlInput_path_name"))
{
  // All generated meshes are regular orthogonal meshes
  _regular_orthogonal_mesh = true;

  // We support 2D or 3D mesh only at this point. If you need 1D
  // Please contact MOOSE developers
  if (_dim == 1)
    mooseError("Support 2 or 3 dimensional mesh only");

  /*
    ======= Xolotl related routines below =======
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

  _xolotl_local_index_table = init_xolotl_local_index_table(6);
  fillout_xolotl_local_index_table(_xolotl_local_index_table, 6);
  //print_xolotl_local_index_table(_xolotl_local_index_table, 6);

  _xolotl_xi_lb = _xolotl_local_index_table[_moose_rank][0];
  _xolotl_xi_ub = _xolotl_xi_lb + _xolotl_local_index_table[_moose_rank][1] - 1;
  _xolotl_yi_lb = _xolotl_local_index_table[_moose_rank][2];
  _xolotl_yi_ub = _xolotl_yi_lb + _xolotl_local_index_table[_moose_rank][3] - 1;
  _xolotl_zi_lb = _xolotl_local_index_table[_moose_rank][4];
  _xolotl_zi_ub = _xolotl_zi_lb + _xolotl_local_index_table[_moose_rank][5] - 1;
  _xolotl_localNx = _xolotl_local_index_table[_moose_rank][1];
  _xolotl_localNy = _xolotl_local_index_table[_moose_rank][3];
  _xolotl_localNz = _xolotl_local_index_table[_moose_rank][5];

  count_num_partition_along_each_xyz(&_xnp, &_ynp, &_znp, _xolotl_local_index_table);
  //printf("_xnp = %d, _ynp = %d, _znp = %d\n", _xnp, _ynp, _znp);

  get_current_process_coord(&_xpid, &_ypid, &_zpid, _moose_rank);

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

  _xolotl_interface->finalizeXolotl(_xolotl_solver, ISSTANDALONE); 

}

void
XolotlMeshSynced::buildMesh()
{
  // Reference to the libmesh mesh
  MeshBase & mesh = getMesh();

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

      default:
        mooseError("Does not support dimension ", _dim, "yet");
    }
  }

  _elem_type = Utility::string_to_enum<ElemType>(elem_type_enum);

  mesh.set_mesh_dimension(_dim);
  mesh.set_spatial_dimension(_dim);

  // Switching on MooseEnum
  switch (_dim)
  {
    case 1:
      build_segment_Edge2(dynamic_cast<UnstructuredMesh &>(getMesh()), _elem_type);
      break;
    case 2:
      build_square_Quad4(dynamic_cast<UnstructuredMesh &>(getMesh()), _elem_type);
      break;
    case 3:
      build_cube_Hex8(dynamic_cast<UnstructuredMesh &>(getMesh()), _elem_type);
      break;
    default:
      mooseError("Does not support dimension ", _dim, "yet");
  }
}

Real
XolotlMeshSynced::getMinInDimension(unsigned int component) const
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
XolotlMeshSynced::getMaxInDimension(unsigned int component) const
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
XolotlMeshSynced::safeClone() const
{
  return libmesh_make_unique<XolotlMeshSynced>(*this);
}

dof_id_type
XolotlMeshSynced::node_id_Edge2(const ElemType /*type*/,
              const dof_id_type /*nx*/,
              const dof_id_type /*ny*/,
              const dof_id_type i,
              const dof_id_type /*j*/,
              const dof_id_type /*k*/) const

{
  // Transform a grid coordinate (i, j) to its global node ID
  // This match what PETSc does
  return i;
}

dof_id_type
XolotlMeshSynced::node_id_Quad4(const ElemType /*type*/,
              const dof_id_type nx,
              const dof_id_type /*ny*/,
              const dof_id_type i,
              const dof_id_type j,
              const dof_id_type /*k*/) const

{
  // Transform a grid coordinate (i, j) to its global node ID
  // This match what PETSc does
  return i + j * (nx + 1);
}

dof_id_type
XolotlMeshSynced::node_id_Hex8(const ElemType /*type*/,
              const dof_id_type nx,
              const dof_id_type ny,
              const dof_id_type nz,
              const dof_id_type i,
              const dof_id_type j,
              const dof_id_type k) const

{
  // Transform a grid coordinate (i, j) to its global node ID
  // This match what PETSc does
  return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
}

void
XolotlMeshSynced::add_element_Edge2(
                  const dof_id_type nx,
                  const dof_id_type i,
                  const dof_id_type elem_id,
                  const processor_id_type pid,
                  const ElemType type,
                  MeshBase & mesh) const
{
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  int xp, xpid, xpidplus;
  xp = _xnp;
  xpidplus = _xpid;
  xpid = (i == _xolotl_xi_lb - 1) ? xpidplus - 1 : xpidplus;

  // Left
  Real xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  auto node0_ptr = mesh.add_point(libMesh::Point(xcoord, 0, 0),
                                  node_id_Edge2(type, 0, 0, i, 0, 0));
  node0_ptr->set_unique_id() = node_id_Edge2(type, 0, 0, i, 0, 0);
  node0_ptr->set_id() = node0_ptr->unique_id();
  // xpid + ypid * xp is the global processor ID
  node0_ptr->processor_id() = xpid;

  // Right
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  auto node1_ptr =
      mesh.add_point(libMesh::Point(xcoord, 0, 0),
                     node_id_Quad4(type, 0, 0, i + 1, 0, 0));
  node1_ptr->set_unique_id() = node_id_Quad4(type, 0, 0, i + 1, 0, 0);
  node1_ptr->set_id() = node1_ptr->unique_id();
  node1_ptr->processor_id() = xpidplus;

  // New an element and attach four nodes to it
  Elem * elem = new Edge2;
  elem->set_id(elem_id);
  elem->processor_id() = pid;
  elem->set_unique_id() = elem_id;
  elem = mesh.add_elem(elem);
  elem->set_node(0) = node0_ptr;
  elem->set_node(1) = node1_ptr;

  // Left
  if (i == 0)
    boundary_info.add_side(elem, 0, 0);

  // Right
  if (i == nx - 1)
    boundary_info.add_side(elem, 1, 1);
}

void
// XolotlMeshSynced::add_element_Quad4(DM da,
XolotlMeshSynced::add_element_Quad4(
                  const dof_id_type nx,
                  const dof_id_type ny,
                  const dof_id_type i,
                  const dof_id_type j,
                  const dof_id_type elem_id,
                  const processor_id_type pid,
                  const ElemType type,
                  MeshBase & mesh) const
{
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  int xp, yp, xpid, ypid, xpidplus, ypidplus;
  xp = _xnp;
  yp = _ynp;
  xpidplus = _xpid;
  xpid = (i == _xolotl_xi_lb - 1) ? xpidplus - 1 : xpidplus;
  ypidplus = _ypid;
  ypid = (j == _xolotl_yi_lb - 1) ? ypidplus - 1 : ypidplus;

  // Bottom Left
  Real xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  Real ycoord = static_cast<Real> (_xolotl_yc[j]);
  auto node0_ptr = mesh.add_point(libMesh::Point(xcoord, ycoord, 0),
                                  node_id_Quad4(type, nx, 0, i, j, 0));
  node0_ptr->set_unique_id() = node_id_Quad4(type, nx, 0, i, j, 0);
  node0_ptr->set_id() = node0_ptr->unique_id();
  // xpid + ypid * xp is the global processor ID
  node0_ptr->processor_id() = xpid + ypid * xp;

  // Bottom Right
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  ycoord = static_cast<Real> (_xolotl_yc[j]);
  auto node1_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, 0),
                     node_id_Quad4(type, nx, 0, i + 1, j, 0));
  node1_ptr->set_unique_id() = node_id_Quad4(type, nx, 0, i + 1, j, 0);
  node1_ptr->set_id() = node1_ptr->unique_id();
  node1_ptr->processor_id() = xpidplus + ypid * xp;

  // Top Right
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  ycoord = static_cast<Real> (_xolotl_yc[j + 1]);
  auto node2_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, 0),
                     node_id_Quad4(type, nx, 0, i + 1, j + 1, 0));
  node2_ptr->set_unique_id() = node_id_Quad4(type, nx, 0, i + 1, j + 1, 0);
  node2_ptr->set_id() = node2_ptr->unique_id();
  node2_ptr->processor_id() = xpidplus + ypidplus * xp;

  // Top Left
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  ycoord = static_cast<Real> (_xolotl_yc[j + 1]);
  auto node3_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, 0),
                     node_id_Quad4(type, nx, 0, i, j + 1, 0));
  node3_ptr->set_unique_id() = node_id_Quad4(type, nx, 0, i, j + 1, 0);
  node3_ptr->set_id() = node3_ptr->unique_id();
  node3_ptr->processor_id() = xpid + ypidplus * xp;

  // New an element and attach four nodes to it
  Elem * elem = new Quad4;
  elem->set_id(elem_id);
  elem->processor_id() = pid;
  elem->set_unique_id() = elem_id;
  elem = mesh.add_elem(elem);
  elem->set_node(0) = node0_ptr;
  elem->set_node(1) = node1_ptr;
  elem->set_node(2) = node2_ptr;
  elem->set_node(3) = node3_ptr;

  // Bottom
  if (j == 0)
    boundary_info.add_side(elem, 0, 0);

  // Right
  if (i == nx - 1)
    boundary_info.add_side(elem, 1, 1);

  // Top
  if (j == ny - 1)
    boundary_info.add_side(elem, 2, 2);

  // Left
  if (i == 0)
    boundary_info.add_side(elem, 3, 3);
}

void
// XolotlMeshSynced::add_element_Hex8(DM da,
  XolotlMeshSynced::add_element_Hex8(
                  const dof_id_type nx,
                  const dof_id_type ny,
                  const dof_id_type nz,
                  const dof_id_type i,
                  const dof_id_type j,
                  const dof_id_type k,
                  const dof_id_type elem_id,
                  const processor_id_type pid,
                  const ElemType type,
                  MeshBase & mesh) const
{
  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  int xp, yp, zp, xpid, ypid, zpid, xpidplus, ypidplus, zpidplus;
  xp = _xnp;
  yp = _ynp;
  zp = _znp;
  xpidplus = _xpid;
  xpid = (i == _xolotl_xi_lb - 1) ? xpidplus - 1 : xpidplus;
  ypidplus = _ypid;
  ypid = (j == _xolotl_yi_lb - 1) ? ypidplus - 1 : ypidplus;
  zpidplus = _zpid;
  zpid = (k == _xolotl_zi_lb - 1) ? zpidplus - 1 : zpidplus;

  // Bottom Left, back
  Real xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  Real ycoord = static_cast<Real> (_xolotl_yc[j]);
  Real zcoord = static_cast<Real> (_xolotl_zc[k]);
  auto node0_ptr = mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                                  node_id_Hex8(type, nx, ny, 0, i, j, k));
  node0_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i, j, k);
  node0_ptr->set_id() = node0_ptr->unique_id();
  // xpid + ypid * xp is the global processor ID
  node0_ptr->processor_id() = xpid + ypid * xp + zpid * xp * yp;

  // Bottom Right, back
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  ycoord = static_cast<Real> (_xolotl_yc[j]);
  zcoord = static_cast<Real> (_xolotl_zc[k]);
  auto node1_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                     node_id_Hex8(type, nx, ny, 0, i + 1, j, k));
  node1_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i + 1, j, k);
  node1_ptr->set_id() = node1_ptr->unique_id();
  node1_ptr->processor_id() = xpidplus + ypid * xp + zpid * xp * yp;

  // Top Right, back
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  ycoord = static_cast<Real> (_xolotl_yc[j + 1]);
  zcoord = static_cast<Real> (_xolotl_zc[k]);
  auto node2_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                     node_id_Hex8(type, nx, ny, 0, i + 1, j + 1, k));
  node2_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i + 1, j + 1, k);
  node2_ptr->set_id() = node2_ptr->unique_id();
  node2_ptr->processor_id() = xpidplus + ypidplus * xp + zpid * xp * yp;

  // Top Left, back
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  ycoord = static_cast<Real> (_xolotl_yc[j + 1]);
  zcoord = static_cast<Real> (_xolotl_zc[k]);
  auto node3_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                     node_id_Hex8(type, nx, ny, 0, i, j + 1, k));
  node3_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i, j + 1, k);
  node3_ptr->set_id() = node3_ptr->unique_id();
  node3_ptr->processor_id() = xpid + ypidplus * xp + zpid * xp * yp;

  // Bottom Left, front
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  ycoord = static_cast<Real> (_xolotl_yc[j]);
  zcoord = static_cast<Real> (_xolotl_zc[k + 1]);
  auto node4_ptr = mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                                  node_id_Hex8(type, nx, ny, 0, i, j, k + 1));
  node4_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i, j, k + 1);
  node4_ptr->set_id() = node4_ptr->unique_id();
  // xpid + ypid * xp is the global processor ID
  node4_ptr->processor_id() = xpid + ypid * xp + zpidplus * xp * yp;

  // Bottom Right, front
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  ycoord = static_cast<Real> (_xolotl_yc[j]);
  zcoord = static_cast<Real> (_xolotl_zc[k + 1]);
  auto node5_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                     node_id_Hex8(type, nx, ny, 0, i + 1, j, k + 1));
  node5_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i + 1, j, k + 1);
  node5_ptr->set_id() = node5_ptr->unique_id();
  node5_ptr->processor_id() = xpidplus + ypid * xp + zpidplus * xp * yp;

  // Top Right, front
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i + 1] : _xolotl_xcNR[i + 1]);
  ycoord = static_cast<Real> (_xolotl_yc[j + 1]);
  zcoord = static_cast<Real> (_xolotl_zc[k + 1]);
  auto node6_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                     node_id_Hex8(type, nx, ny, 0, i + 1, j + 1, k + 1));
  node6_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i + 1, j + 1, k + 1);
  node6_ptr->set_id() = node6_ptr->unique_id();
  node6_ptr->processor_id() = xpidplus + ypidplus * xp + zpidplus * xp * yp;

  // Top Left, front
  xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  ycoord = static_cast<Real> (_xolotl_yc[j + 1]);
  zcoord = static_cast<Real> (_xolotl_zc[k + 1]);
  auto node7_ptr =
      mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                     node_id_Hex8(type, nx, ny, 0, i, j + 1, k + 1));
  node7_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, 0, i, j + 1, k + 1);
  node7_ptr->set_id() = node7_ptr->unique_id();
  node7_ptr->processor_id() = xpid + ypidplus * xp + zpidplus * xp * yp;

  // New an element and attach four nodes to it
  Elem * elem = new Hex8;
  elem->set_id(elem_id);
  elem->processor_id() = pid;
  elem->set_unique_id() = elem_id;
  elem = mesh.add_elem(elem);
  elem->set_node(0) = node0_ptr;
  elem->set_node(1) = node1_ptr;
  elem->set_node(2) = node2_ptr;
  elem->set_node(3) = node3_ptr;
  elem->set_node(4) = node4_ptr;
  elem->set_node(5) = node5_ptr;
  elem->set_node(6) = node6_ptr;
  elem->set_node(7) = node7_ptr;

  // Bottom
  if (j == 0)
    boundary_info.add_side(elem, 0, 0);

  // Right
  if (i == nx - 1)
    boundary_info.add_side(elem, 1, 1);

  // Top
  if (j == ny - 1)
    boundary_info.add_side(elem, 2, 2);

  // Left
  if (i == 0)
    boundary_info.add_side(elem, 3, 3);

  // Front
  if (k == nz - 1)
    boundary_info.add_side(elem, 4, 4);

  // Back
  if (k == 0)
    boundary_info.add_side(elem, 5, 5);
}

void
XolotlMeshSynced::set_boundary_names_Edge2(BoundaryInfo & boundary_info) const
{
  boundary_info.sideset_name(0) = "left";
  boundary_info.sideset_name(1) = "right";
}

void
XolotlMeshSynced::set_boundary_names_Quad4(BoundaryInfo & boundary_info) const
{
  boundary_info.sideset_name(0) = "bottom";
  boundary_info.sideset_name(1) = "right";
  boundary_info.sideset_name(2) = "top";
  boundary_info.sideset_name(3) = "left";
}

void
XolotlMeshSynced::set_boundary_names_Hex8(BoundaryInfo & boundary_info) const
{
  boundary_info.sideset_name(0) = "bottom";
  boundary_info.sideset_name(1) = "right";
  boundary_info.sideset_name(2) = "top";
  boundary_info.sideset_name(3) = "left";
  boundary_info.sideset_name(4) = "front";
  boundary_info.sideset_name(5) = "back";
}

void
XolotlMeshSynced::get_indices_Edge2(const dof_id_type nx,
                  const dof_id_type /*ny*/,
                  const dof_id_type elem_id,
                  dof_id_type & i) const
{
  i = elem_id % nx;
}

void
XolotlMeshSynced::get_indices_Quad4(const dof_id_type nx,
                  const dof_id_type /*ny*/,
                  const dof_id_type elem_id,
                  dof_id_type & i,
                  dof_id_type & j) const
{
  i = elem_id % nx;
  j = (elem_id - i) / nx;
}

void
XolotlMeshSynced::get_indices_Hex8(const dof_id_type nx,
                  const dof_id_type ny,
                  const dof_id_type elem_id,
                  dof_id_type & i,
                  dof_id_type & j,
                  dof_id_type & k) const
{
  i = elem_id % nx;
  j = ((elem_id - i) / nx) % ny;
  k = elem_id / (nx * ny);
}

dof_id_type
XolotlMeshSynced::elem_id_Edge2(const dof_id_type i) const
{
  return i;
}

dof_id_type
XolotlMeshSynced::elem_id_Quad4(const dof_id_type nx,
              const dof_id_type /*nx*/,
              const dof_id_type i,
              const dof_id_type j,
              const dof_id_type /*k*/) const
{
  return (j * nx) + i;
}

dof_id_type
XolotlMeshSynced::elem_id_Hex8(const dof_id_type nx,
              const dof_id_type ny,
              const dof_id_type i,
              const dof_id_type j,
              const dof_id_type k) const
{
  return (k * nx * ny) + (j * nx) + i;
}

void
XolotlMeshSynced::get_neighbors_Edge2(const dof_id_type nx,
                    const dof_id_type i,
                    std::vector<dof_id_type> & neighbors) const
{
  std::fill(neighbors.begin(), neighbors.end(), Elem::invalid_id);


  // Right
  if (i != nx - 1)
    neighbors[1] = elem_id_Edge2(i + 1);

  // Left
  if (i != 0)
    neighbors[0] = elem_id_Edge2(i - 1);
}

void
XolotlMeshSynced::get_neighbors_Quad4(const dof_id_type nx,
                    const dof_id_type ny,
                    const dof_id_type i,
                    const dof_id_type j,
                    std::vector<dof_id_type> & neighbors) const
{
  std::fill(neighbors.begin(), neighbors.end(), Elem::invalid_id);

  // Bottom
  if (j != 0)
    neighbors[0] = elem_id_Quad4(nx, 0, i, j - 1, 0);

  // Right
  if (i != nx - 1)
    neighbors[1] = elem_id_Quad4(nx, 0, i + 1, j, 0);

  // Top
  if (j != ny - 1)
    neighbors[2] = elem_id_Quad4(nx, 0, i, j + 1, 0);

  // Left
  if (i != 0)
    neighbors[3] = elem_id_Quad4(nx, 0, i - 1, j, 0);
}

void
XolotlMeshSynced::get_neighbors_Hex8(const dof_id_type nx,
                    const dof_id_type ny,
                    const dof_id_type nz,
                    const dof_id_type i,
                    const dof_id_type j,
                    const dof_id_type k,
                    std::vector<dof_id_type> & neighbors) const
{
  std::fill(neighbors.begin(), neighbors.end(), Elem::invalid_id);

  // Bottom
  if (j != 0)
    neighbors[0] = elem_id_Hex8(nx, ny, i, j - 1, k);

  // Right
  if (i != nx - 1)
    neighbors[1] = elem_id_Hex8(nx, ny, i + 1, j, k);

  // Top
  if (j != ny - 1)
    neighbors[2] = elem_id_Hex8(nx, ny, i, j + 1, k);

  // Left
  if (i != 0)
    neighbors[3] = elem_id_Hex8(nx, ny, i - 1, j, k);

  // front
  if (k != nz - 1)
    neighbors[4] = elem_id_Hex8(nx, ny, i, j, k + 1);

  // back
  if (k != 0)
    neighbors[5] = elem_id_Hex8(nx, ny, i, j, k - 1);
}

void
XolotlMeshSynced::get_ghost_neighbors_Edge2(const dof_id_type nx,
                          const MeshBase & mesh,
                          std::set<dof_id_type> & ghost_elems) const
{
  auto & boundary_info = mesh.get_boundary_info();

  dof_id_type i;

  std::vector<dof_id_type> neighbors(2);

  for (auto elem_ptr : mesh.element_ptr_range())
  {
    for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
    {
      // No current neighbor
      if (!elem_ptr->neighbor_ptr(s))
      {
        // Not on a boundary
        if (!boundary_info.n_boundary_ids(elem_ptr, s))
        {
          auto elem_id = elem_ptr->id();

          get_indices_Edge2(nx, 0, elem_id, i);

          get_neighbors_Edge2(nx, i, neighbors);

          ghost_elems.insert(neighbors[s]);
        }
      }
    }
  }
}

void
XolotlMeshSynced::get_ghost_neighbors_Quad4(const dof_id_type nx,
                          const dof_id_type ny,
                          const MeshBase & mesh,
                          std::set<dof_id_type> & ghost_elems) const
{
  auto & boundary_info = mesh.get_boundary_info();

  dof_id_type i, j;

  std::vector<dof_id_type> neighbors(4);

  for (auto elem_ptr : mesh.element_ptr_range())
  {
    for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
    {
      // No current neighbor
      if (!elem_ptr->neighbor_ptr(s))
      {
        // Not on a boundary
        if (!boundary_info.n_boundary_ids(elem_ptr, s))
        {
          auto elem_id = elem_ptr->id();

          get_indices_Quad4(nx, 0, elem_id, i, j);

          get_neighbors_Quad4(nx, ny, i, j, neighbors);

          ghost_elems.insert(neighbors[s]);
        }
      }
    }
  }
}

void
XolotlMeshSynced::get_ghost_neighbors_Hex8(const dof_id_type nx,
                          const dof_id_type ny,
                          const dof_id_type nz,
                          const MeshBase & mesh,
                          std::set<dof_id_type> & ghost_elems)  const
{
  auto & boundary_info = mesh.get_boundary_info();

  dof_id_type i, j, k;

  std::vector<dof_id_type> neighbors(6);

  for (auto elem_ptr : mesh.element_ptr_range())
  {
    for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
    {
      // No current neighbor
      if (!elem_ptr->neighbor_ptr(s))
      {
        // Not on a boundary
        if (!boundary_info.n_boundary_ids(elem_ptr, s))
        {
          auto elem_id = elem_ptr->id();

          get_indices_Hex8(nx, ny, elem_id, i, j, k);

          get_neighbors_Hex8(nx, ny, nz, i, j, k, neighbors);

          ghost_elems.insert(neighbors[s]);
        }
      }
    }
  }
}

void
XolotlMeshSynced::add_node_Edg2(dof_id_type nx,
              dof_id_type i,
              processor_id_type pid,
              ElemType type,
              MeshBase & mesh) const
{

  Real xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  auto node0_ptr = mesh.add_point(libMesh::Point(xcoord, 0, 0),
                                  node_id_Edge2(type, 0, 0, i, 0, 0));
  node0_ptr->set_unique_id() = node_id_Edge2(type, 0, 0, i, 0, 0);
  node0_ptr->processor_id() = pid;
}

void
XolotlMeshSynced::add_node_Qua4(dof_id_type nx,
              dof_id_type ny,
              dof_id_type i,
              dof_id_type j,
              processor_id_type pid,
              ElemType type,
              MeshBase & mesh) const
{

  Real xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  Real ycoord = static_cast<Real> (_xolotl_yc[j]);
  auto node0_ptr = mesh.add_point(libMesh::Point(xcoord, ycoord, 0),
                                  node_id_Quad4(type, nx, 0, i, j, 0));
  node0_ptr->set_unique_id() = node_id_Quad4(type, nx, 0, i, j, 0);
  node0_ptr->processor_id() = pid;
}

void
XolotlMeshSynced::add_node_Hex8(dof_id_type nx,
              dof_id_type ny,
              dof_id_type nz,
              dof_id_type i,
              dof_id_type j,
              dof_id_type k,
              processor_id_type pid,
              ElemType type,
              MeshBase & mesh) const
{

  // Bottom Left
  Real xcoord = static_cast<Real> (_xolotl_regulargrid ? _xolotl_xc[i] : _xolotl_xcNR[i]);
  Real ycoord = static_cast<Real> (_xolotl_yc[j]);
  Real zcoord = static_cast<Real> (_xolotl_zc[k]);
  auto node0_ptr = mesh.add_point(libMesh::Point(xcoord, ycoord, zcoord),
                                  node_id_Hex8(type, nx, ny, nz, i, j, k));
  node0_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, nz, i, j, k);
  node0_ptr->processor_id() = pid;
}

void
XolotlMeshSynced::build_segment_Edge2(UnstructuredMesh & mesh, const ElemType type) const
{
  const auto pid = mesh.comm().rank();

  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  int xs, xm, Mx;
  xs = _xolotl_xi_lb;
  xm = _xolotl_localNx;
  Mx = _xolotl_nx;

  for (PetscInt i = xs; i < xs + xm; i++)
  {
    // We loop over grid points, but we are
    // here building elements. So that we just
    // simply skip the first x points since the
    // number of grid ponts is one more than
    // the number of grid elements
    if (!i)
      continue;

    dof_id_type ele_id = (i - 1);

    add_element_Edge2(Mx - 1, i - 1, ele_id, pid, type, mesh);
  }

  // If there is no element at the given processor
  // We need to manually add all mesh nodes
  if ((xs == 0 && xm == 1))
    for (PetscInt i = xs; i < xs + xm; i++)
      add_node_Edg2(Mx, i, pid, type, mesh);

  // Need to link up the local elements before we can know what's missing
  mesh.find_neighbors();

  mesh.find_neighbors(true);

  // Set RemoteElem neighbors
  for (auto & elem_ptr : mesh.element_ptr_range())
    for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
      if (!elem_ptr->neighbor_ptr(s) && !boundary_info.n_boundary_ids(elem_ptr, s))
        elem_ptr->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

  set_boundary_names_Quad4(boundary_info);
  // Already partitioned!
  mesh.skip_partitioning(true);

  // No need to renumber or find neighbors - done did it.
  // Avoid deprecation message/error by _also_ setting
  // allow_renumbering(false). This is a bit silly, but we want to
  // catch cases where people are purely using the old "skip"
  // interface and not the new flag setting one.
  mesh.allow_renumbering(false);
  mesh.prepare_for_use(/*skip_renumber (ignored!) = */ false,
                       /*skip_find_neighbors = */ true);
}

void
XolotlMeshSynced::build_square_Quad4(UnstructuredMesh & mesh, const ElemType type) const
{
  const auto pid = mesh.comm().rank();

  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  int xs, ys, xm, ym, Mx, My;
  xs = _xolotl_xi_lb;
  xm = _xolotl_localNx;
  ys = _xolotl_yi_lb;
  ym = _xolotl_localNy;
  Mx = _xolotl_nx;
  My = _xolotl_ny;

  //printf("pid = %d, xs = %d, xm = %d, ys = %d,  ym = %d\n", pid, xs, xm, ys, ym);

  for (PetscInt j = ys; j < ys + ym; j++)
    for (PetscInt i = xs; i < xs + xm; i++)
    {
      // We loop over grid points, but we are
      // here building elements. So that we just
      // simply skip the first x and y points since the
      // number of grid ponts is one more than
      // the number of grid elements
      if (!i || !j)
        continue;

      dof_id_type ele_id = (i - 1) + (j - 1) * (Mx - 1);

      // add_element_Quad4(da, Mx - 1, My - 1, i - 1, j - 1, ele_id, pid, type, mesh);
      add_element_Quad4(Mx - 1, My - 1, i - 1, j - 1, ele_id, pid, type, mesh);
    }

  // If there is no element at the given processor
  // We need to manually add all mesh nodes
  if ((ys == 0 && ym == 1) || (xs == 0 && xm == 1))
    for (PetscInt j = ys; j < ys + ym; j++)
      for (PetscInt i = xs; i < xs + xm; i++)
        add_node_Qua4(Mx, My, i, j, pid, type, mesh);

  // Need to link up the local elements before we can know what's missing
  mesh.find_neighbors();

  mesh.find_neighbors(true);

  // Set RemoteElem neighbors
  for (auto & elem_ptr : mesh.element_ptr_range())
    for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
      if (!elem_ptr->neighbor_ptr(s) && !boundary_info.n_boundary_ids(elem_ptr, s))
        elem_ptr->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

  set_boundary_names_Quad4(boundary_info);
  // Already partitioned!
  mesh.skip_partitioning(true);

  // No need to renumber or find neighbors - done did it.
  // Avoid deprecation message/error by _also_ setting
  // allow_renumbering(false). This is a bit silly, but we want to
  // catch cases where people are purely using the old "skip"
  // interface and not the new flag setting one.
  mesh.allow_renumbering(false);
  mesh.prepare_for_use(/*skip_renumber (ignored!) = */ false,
                       /*skip_find_neighbors = */ true);
}

void
XolotlMeshSynced::build_cube_Hex8(UnstructuredMesh & mesh, const ElemType type) const
{
  const auto pid = mesh.comm().rank();

  BoundaryInfo & boundary_info = mesh.get_boundary_info();

  int xs, xm, ys, ym, zs, zm, Mx, My, Mz;
  xs = _xolotl_xi_lb;
  xm = _xolotl_localNx;
  ys = _xolotl_yi_lb;
  ym = _xolotl_localNy;
  zs = _xolotl_zi_lb;
  zm = _xolotl_localNz;
  Mx = _xolotl_nx;
  My = _xolotl_ny;
  Mz = _xolotl_nz;

  for (PetscInt k = zs; k < zs + zm; k++)
    for (PetscInt j = ys; j < ys + ym; j++)
      for (PetscInt i = xs; i < xs + xm; i++)
    {
      // We loop over grid points, but we are
      // here building elements. So that we just
      // simply skip the first x and y points since the
      // number of grid ponts is one more than
      // the number of grid elements
      if (!i || !j || !k)
        continue;

      dof_id_type ele_id = (i - 1) + (j - 1) * (Mx - 1) + (k - 1) * (Mx - 1) * (My - 1);

      add_element_Hex8(Mx - 1, My - 1, Mz - 1, i - 1, j - 1, k - 1, ele_id, pid, type, mesh);
    }

    // done here

  // If there is no element at the given processor
  // We need to manually add all mesh nodes
  if ((ys == 0 && ym == 1) || (xs == 0 && xm == 1) || (zs == 0 && zm == 1))
    for (PetscInt k = zs; k < zs + zm; k++)
      for (PetscInt j = ys; j < ys + ym; j++)
        for (PetscInt i = xs; i < xs + xm; i++)
          // add_node_Qua4(Mx, My, i, j, pid, type, mesh, interface);
          add_node_Hex8(Mx, My, Mz, i, j, k, pid, type, mesh);

  // Need to link up the local elements before we can know what's missing
  mesh.find_neighbors();

  mesh.find_neighbors(true);

  // Set RemoteElem neighbors
  for (auto & elem_ptr : mesh.element_ptr_range())
    for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
      if (!elem_ptr->neighbor_ptr(s) && !boundary_info.n_boundary_ids(elem_ptr, s))
        elem_ptr->set_neighbor(s, const_cast<RemoteElem *>(remote_elem));

  set_boundary_names_Hex8(boundary_info);
  // Already partitioned!
  mesh.skip_partitioning(true);

  // No need to renumber or find neighbors - done did it.
  // Avoid deprecation message/error by _also_ setting
  // allow_renumbering(false). This is a bit silly, but we want to
  // catch cases where people are purely using the old "skip"
  // interface and not the new flag setting one.
  mesh.allow_renumbering(false);
  mesh.prepare_for_use(/*skip_renumber (ignored!) = */ false,
                       /*skip_find_neighbors = */ true);
}

std::vector<double>
XolotlMeshSynced::build_xolotl_axis(int nsize, double dl) const
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

int**
XolotlMeshSynced::init_xolotl_local_index_table(int ncolumn) const
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
XolotlMeshSynced::print_xolotl_local_index_table(int **table, int ncol) const
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
XolotlMeshSynced::fillout_xolotl_local_index_table(int **table, int ncol) const
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

void
XolotlMeshSynced::count_num_partition_along_each_xyz(int *xp, int *yp, int *zp, int **local_index_table) const
{
  int nrow;
  MPI_Comm_size(MPI_COMM_WORLD, &nrow);
  int numx = 0;
  int numy = 0;
  int numz = 0;
  for (int i = 0; i < nrow; i++) {
    int xs = local_index_table[i][0];
    int ys = local_index_table[i][2];
    int zs = local_index_table[i][4];
    //Searching through x axis
    if (ys == 0 && zs == 0) {
      numx++;
    }
    //Searching through y axis
    if (zs == 0 && xs == 0) {
      numy++;
    }
    //Searching through z axis
    if (xs == 0 && ys == 0) {
      numz++;
    }
  }
  *xp = numx;
  *yp = numy;
  *zp = numz;
}

void
XolotlMeshSynced::get_current_process_coord(int *xpid, int *ypid, int *zpid, int pid) const
{
  // pid = xpid + ypid * _xnp + zpid * _xnp * _ynp
  *xpid = pid % _xnp;
  *ypid = (pid / _xnp) % _ynp;
  *zpid = pid / (_xnp * _ynp);
}
