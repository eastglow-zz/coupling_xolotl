//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PETScDMDAMesh.h"

#include "libmesh/mesh_generation.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/partitioner.h"
#include "libmesh/metis_csr_graph.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/remote_elem.h"
#include "libmesh/face_quad4.h"
#include "libmesh/cell_hex8.h"

#include "coupling_xolotlApp.h"
// C++ includes
#include <cmath> // provides round, not std::round (see http://www.cplusplus.com/reference/cmath/round/)

registerMooseObject("coupling_xolotlApp", PETScDMDAMesh);

template<>
InputParameters validParams<PETScDMDAMesh>() {
	InputParameters params = validParams<MooseMesh>();

	MooseEnum elem_types("EDGE2  QUAD4  HEX8"); // no default

	MooseEnum dims("1=1 2 3", "2");
	params.addRequiredParam < MooseEnum
			> ("dim", dims, "The dimension of the mesh to be generated"); // Make this parameter required

	params.addParam < MooseEnum
			> ("elem_type", elem_types, "The type of element from libMesh to "
					"generate (default: linear element for "
					"requested dimension)");

	params.addParamNamesToGroup("dim", "Main");

	params.addClassDescription(
			"Create a line, square, or cube mesh with uniformly spaced memsh using PETSc DMDA.");

	// This mesh is always distributed
	params.set < MooseEnum > ("parallel_type") = "DISTRIBUTED";

	// Parameter for the Xolotl file name
	params.addRequiredParam < FileName
			> ("XolotlInput_path_name", "Name with the path for the Xolotl input file");

	return params;
}

PETScDMDAMesh::PETScDMDAMesh(const InputParameters & parameters) :
		MooseMesh(parameters), _xolotl_input_path_name(
				getParam < FileName > ("XolotlInput_path_name")), _dim(
				getParam < MooseEnum > ("dim")) {
	// All generated meshes are regular orthogonal meshes
	_regular_orthogonal_mesh = true;

	// Get the external app to create the interface and its grid
	coupling_xolotlApp * xolotl_app = dynamic_cast<coupling_xolotlApp *>(&_app);
	if (xolotl_app) {
		auto & interface = xolotl_app->getInterface();
		// This has to be done here because the base app cannot take parameters
		// from the input file
		int argc = 2;
		char ** argv = new char*[argc];
		std::string parameterFile = "bla";
		argv[0] = new char[parameterFile.length() + 1];
		strcpy(argv[0], parameterFile.c_str());
		argv[1] = new char[_xolotl_input_path_name.length() + 1];
		strcpy(argv[1], _xolotl_input_path_name.c_str());

		interface.initializeXolotl(argc, argv, _communicator.get(), false);
		// Now we can get the TS from the app
		TS & ts = xolotl_app->getXolotlTS();
		// Retrieve mesh from TS
		TSGetDM(ts, &_dmda);
	} else
		mooseError("Missing the Xolotl App");
}

Real PETScDMDAMesh::getMinInDimension(unsigned int component) const {
	switch (component) {
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

Real PETScDMDAMesh::getMaxInDimension(unsigned int component) const {
	switch (component) {
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

std::unique_ptr<MooseMesh> PETScDMDAMesh::safeClone() const {
	return libmesh_make_unique < PETScDMDAMesh > (*this);
}

inline dof_id_type node_id_Edge2(const ElemType /*type*/, const dof_id_type i)

{
	// Transform a grid coordinate (i, j, k) to its global node ID
	// This match what PETSc does
	return i;
}

inline dof_id_type node_id_Quad4(const ElemType /*type*/, const dof_id_type nx,
		const dof_id_type i, const dof_id_type j)

		{
	// Transform a grid coordinate (i, j, k) to its global node ID
	// This match what PETSc does
	return i + j * (nx + 1);
}

inline dof_id_type node_id_Hex8(const ElemType /*type*/, const dof_id_type nx,
		const dof_id_type ny, const dof_id_type i, const dof_id_type j,
		const dof_id_type k)

		{
	// Transform a grid coordinate (i, j, k) to its global node ID
	// This match what PETSc does
	return i + (j + (ny + 1) * k) * (nx + 1);
}

void add_element_Edge2(DM da, const dof_id_type nx, const dof_id_type i,
		const dof_id_type elem_id, const processor_id_type pid,
		const ElemType type, MeshBase & mesh, XolotlInterface & interface) {
	BoundaryInfo & boundary_info = mesh.get_boundary_info();
	// Mx: number of grid points in x direction for all processors
	// xp: number of processors in x direction
	PetscInt Mx, xp;
	DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE, &xp,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);

	const PetscInt *lx;
	PetscInt *lxo;
	DMGetWorkArray(da, xp + 2, MPIU_INT, &lxo);
	// Gets the ranges of indices in the x, y and z direction that are owned by each process
	// Ranges here are different from what we have in Mat and Vec.
	// It means how many points each processor holds
	DMDAGetOwnershipRanges(da, &lx, PETSC_IGNORE, PETSC_IGNORE);
	lxo[0] = 0;
	for (PetscInt m = 0; m < xp; m++)
		lxo[m + 1] = lxo[m] + lx[m];

	// Try to calculate processor-grid coordinate (xpid, ypid, zpid)
	PetscInt xpid, xpidplus;
	// Finds integer in a sorted array of integers
	// Loc:  the location if found, otherwise -(slot+1)
	// where slot is the place the value would go
	PetscFindInt(i, xp + 1, lxo, &xpid);
	xpid = xpid < 0 ? -xpid - 1 - 1 : xpid;
	PetscFindInt(i + 1, xp + 1, lxo, &xpidplus);
	xpidplus = xpidplus < 0 ? -xpidplus - 1 - 1 : xpidplus;

	DMRestoreWorkArray(da, xp + 2, MPIU_INT, &lxo);

	// Get the geometry of the Xolotl grid information
	double hy = 0.0, hz = 0.0;
	auto xolotlGrid = interface.getGridInfo(hy, hz);

	// Left
	auto node0_ptr = mesh.add_point(libMesh::Point(xolotlGrid[i], 0, 0),
			node_id_Edge2(type, i));
	node0_ptr->set_unique_id() = node_id_Edge2(type, i);
	node0_ptr->set_id() = node0_ptr->unique_id();
	// xpid + ypid * xp is the global processor ID
	node0_ptr->processor_id() = xpid;

	// Right
	auto node1_ptr = mesh.add_point(libMesh::Point(xolotlGrid[i + 1], 0, 0),
			node_id_Edge2(type, i + 1));
	node1_ptr->set_unique_id() = node_id_Edge2(type, i + 1);
	node1_ptr->set_id() = node1_ptr->unique_id();
	node1_ptr->processor_id() = xpidplus;

	// New an element and attach two nodes to it
	Elem * elem = new Edge2;
	elem->set_id(elem_id);
	elem->processor_id() = pid;
	elem->set_unique_id() = elem_id;
	elem = mesh.add_elem(elem);
	elem->set_node(0) = node0_ptr;
	elem->set_node(1) = node1_ptr;

	// Right
	if (i == nx - 1)
		boundary_info.add_side(elem, 0, 0);
	// Left
	if (i == 0)
		boundary_info.add_side(elem, 1, 1);
}

void add_element_Quad4(DM da, const dof_id_type nx, const dof_id_type ny,
		const dof_id_type i, const dof_id_type j, const dof_id_type elem_id,
		const processor_id_type pid, const ElemType type, MeshBase & mesh,
		XolotlInterface & interface) {
	BoundaryInfo & boundary_info = mesh.get_boundary_info();
	// Mx: number of grid points in x direction for all processors
	// My: number of grid points in y direction for all processors
	// xp: number of processors in x direction
	// yp: number of processors in y direction
	PetscInt Mx, My, xp, yp;
	DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, &xp, &yp,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);

	const PetscInt *lx, *ly;
	PetscInt *lxo, *lyo;
	DMGetWorkArray(da, xp + yp + 2, MPIU_INT, &lxo);
	// Gets the ranges of indices in the x, y and z direction that are owned by each process
	// Ranges here are different from what we have in Mat and Vec.
	// It means how many points each processor holds
	DMDAGetOwnershipRanges(da, &lx, &ly, PETSC_IGNORE);
	lxo[0] = 0;
	for (PetscInt m = 0; m < xp; m++)
		lxo[m + 1] = lxo[m] + lx[m];

	lyo = lxo + xp + 1;
	lyo[0] = 0;
	for (PetscInt m = 0; m < yp; m++)
		lyo[m + 1] = lyo[m] + ly[m];

	// Try to calculate processor-grid coordinate (xpid, ypid, zpid)
	PetscInt xpid, ypid, xpidplus, ypidplus;
	// Finds integer in a sorted array of integers
	// Loc:  the location if found, otherwise -(slot+1)
	// where slot is the place the value would go
	PetscFindInt(i, xp + 1, lxo, &xpid);
	xpid = xpid < 0 ? -xpid - 1 - 1 : xpid;
	PetscFindInt(i + 1, xp + 1, lxo, &xpidplus);
	xpidplus = xpidplus < 0 ? -xpidplus - 1 - 1 : xpidplus;

	PetscFindInt(j, yp + 1, lyo, &ypid);
	ypid = ypid < 0 ? -ypid - 1 - 1 : ypid;
	PetscFindInt(j + 1, yp + 1, lyo, &ypidplus);
	ypidplus = ypidplus < 0 ? -ypidplus - 1 - 1 : ypidplus;

	DMRestoreWorkArray(da, xp + yp + 2, MPIU_INT, &lxo);

	// Get the geometry of the Xolotl grid information
	double hy = 0.0, hz = 0.0;
	auto xolotlGrid = interface.getGridInfo(hy, hz);

	// Bottom Left
	auto node0_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j) * hy, 0),
			node_id_Quad4(type, nx, i, j));
	node0_ptr->set_unique_id() = node_id_Quad4(type, nx, i, j);
	node0_ptr->set_id() = node0_ptr->unique_id();
	// xpid + ypid * xp is the global processor ID
	node0_ptr->processor_id() = xpid + ypid * xp;

	// Bottom Right
	auto node1_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i + 1], static_cast<Real>(j) * hy, 0),
			node_id_Quad4(type, nx, i + 1, j));
	node1_ptr->set_unique_id() = node_id_Quad4(type, nx, i + 1, j);
	node1_ptr->set_id() = node1_ptr->unique_id();
	node1_ptr->processor_id() = xpidplus + ypid * xp;

	// Top Right
	auto node2_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i + 1], static_cast<Real>(j + 1) * hy, 0),
			node_id_Quad4(type, nx, i + 1, j + 1));
	node2_ptr->set_unique_id() = node_id_Quad4(type, nx, i + 1, j + 1);
	node2_ptr->set_id() = node2_ptr->unique_id();
	node2_ptr->processor_id() = xpidplus + ypidplus * xp;

	// Top Left
	auto node3_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j + 1) * hy, 0),
			node_id_Quad4(type, nx, i, j + 1));
	node3_ptr->set_unique_id() = node_id_Quad4(type, nx, i, j + 1);
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

void add_element_Hex8(DM da, const dof_id_type nx, const dof_id_type ny,
		const dof_id_type nz, const dof_id_type i, const dof_id_type j,
		const dof_id_type k, const dof_id_type elem_id,
		const processor_id_type pid, const ElemType type, MeshBase & mesh,
		XolotlInterface & interface) {
	BoundaryInfo & boundary_info = mesh.get_boundary_info();
// Mx: number of grid points in x direction for all processors
// My: number of grid points in y direction for all processors
// Mz: number of grid points in z direction for all processors
// xp: number of processors in x direction
// yp: number of processors in y direction
// zp: number of processors in z direction
	PetscInt Mx, My, Mz, xp, yp, zp;
	DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, &xp, &yp, &zp, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);

	const PetscInt *lx, *ly, *lz;
	PetscInt *lxo, *lyo, *lzo;
	DMGetWorkArray(da, xp + yp + zp + 2, MPIU_INT, &lxo);
// Gets the ranges of indices in the x, y and z direction that are owned by each process
// Ranges here are different from what we have in Mat and Vec.
// It means how many points each processor holds
	DMDAGetOwnershipRanges(da, &lx, &ly, &lz);
	lxo[0] = 0;
	for (PetscInt m = 0; m < xp; m++)
		lxo[m + 1] = lxo[m] + lx[m];

	lyo = lxo + xp + 1;
	lyo[0] = 0;
	for (PetscInt m = 0; m < yp; m++)
		lyo[m + 1] = lyo[m] + ly[m];

	lzo = lyo + yp + 1;
	lzo[0] = 0;
	for (PetscInt m = 0; m < zp; m++)
		lzo[m + 1] = lzo[m] + lz[m];

// Try to calculate processor-grid coordinate (xpid, ypid, zpid)
	PetscInt xpid, ypid, zpid, xpidplus, ypidplus, zpidplus;
// Finds integer in a sorted array of integers
// Loc:  the location if found, otherwise -(slot+1)
// where slot is the place the value would go
	PetscFindInt(i, xp + 1, lxo, &xpid);
	xpid = xpid < 0 ? -xpid - 1 - 1 : xpid;
	PetscFindInt(i + 1, xp + 1, lxo, &xpidplus);
	xpidplus = xpidplus < 0 ? -xpidplus - 1 - 1 : xpidplus;

	PetscFindInt(j, yp + 1, lyo, &ypid);
	ypid = ypid < 0 ? -ypid - 1 - 1 : ypid;
	PetscFindInt(j + 1, yp + 1, lyo, &ypidplus);
	ypidplus = ypidplus < 0 ? -ypidplus - 1 - 1 : ypidplus;

	PetscFindInt(k, zp + 1, lzo, &zpid);
	zpid = zpid < 0 ? -zpid - 1 - 1 : zpid;
	PetscFindInt(k + 1, zp + 1, lzo, &zpidplus);
	zpidplus = zpidplus < 0 ? -zpidplus - 1 - 1 : zpidplus;

	DMRestoreWorkArray(da, xp + yp + zp + 2, MPIU_INT, &lxo);

// Get the geometry of the Xolotl grid information
	double hy = 0.0, hz = 0.0;
	auto xolotlGrid = interface.getGridInfo(hy, hz);

// Bottom Left Back
	auto node0_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j) * hy,
					static_cast<Real>(k) * hz),
			node_id_Hex8(type, nx, ny, i, j, k));
	node0_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i, j, k);
	node0_ptr->set_id() = node0_ptr->unique_id();
	node0_ptr->processor_id() = xpid + (ypid + zpid * yp) * xp;

// Bottom Right Back
	auto node1_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i + 1], static_cast<Real>(j) * hy,
					static_cast<Real>(k) * hz),
			node_id_Hex8(type, nx, ny, i + 1, j, k));
	node1_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i + 1, j, k);
	node1_ptr->set_id() = node1_ptr->unique_id();
	node1_ptr->processor_id() = xpidplus + (ypid + zpid * yp) * xp;

// Top Right Back
	auto node2_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i + 1], static_cast<Real>(j + 1) * hy,
					static_cast<Real>(k) * hz),
			node_id_Hex8(type, nx, ny, i + 1, j + 1, k));
	node2_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i + 1, j + 1, k);
	node2_ptr->set_id() = node2_ptr->unique_id();
	node2_ptr->processor_id() = xpidplus + (ypidplus + zpid * yp) * xp;

// Top Left Back
	auto node3_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j + 1) * hy,
					static_cast<Real>(k) * hz),
			node_id_Hex8(type, nx, ny, i, j + 1, k));
	node3_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i, j + 1, k);
	node3_ptr->set_id() = node3_ptr->unique_id();
	node3_ptr->processor_id() = xpid + (ypidplus + zpid * yp) * xp;

// Bottom Left Front
	auto node4_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j) * hy,
					static_cast<Real>(k + 1) * hz),
			node_id_Hex8(type, nx, ny, i, j, k + 1));
	node4_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i, j, k + 1);
	node4_ptr->set_id() = node0_ptr->unique_id();
	node4_ptr->processor_id() = xpid + (ypid + zpidplus * yp) * xp;

// Bottom Right Front
	auto node5_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i + 1], static_cast<Real>(j) * hy,
					static_cast<Real>(k + 1) * hz),
			node_id_Hex8(type, nx, ny, i + 1, j, k + 1));
	node5_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i + 1, j, k + 1);
	node5_ptr->set_id() = node1_ptr->unique_id();
	node5_ptr->processor_id() = xpidplus + (ypid + zpidplus * yp) * xp;

// Top Right Front
	auto node6_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i + 1], static_cast<Real>(j + 1) * hy,
					static_cast<Real>(k + 1) * hz),
			node_id_Hex8(type, nx, ny, i + 1, j + 1, k + 1));
	node6_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i + 1, j + 1,
			k + 1);
	node6_ptr->set_id() = node2_ptr->unique_id();
	node6_ptr->processor_id() = xpidplus + (ypidplus + zpidplus * yp) * xp;

// Top Left Front
	auto node7_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j + 1) * hy,
					static_cast<Real>(k + 1) * hz),
			node_id_Hex8(type, nx, ny, i, j + 1, k + 1));
	node7_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i, j + 1, k + 1);
	node7_ptr->set_id() = node3_ptr->unique_id();
	node7_ptr->processor_id() = xpid + (ypidplus + zpidplus * yp) * xp;

// New an element and attach eight nodes to it
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

// Back
	if (k == 0)
		boundary_info.add_side(elem, 0, 0);
// Bottom
	if (j == 0)
		boundary_info.add_side(elem, 1, 1);
// Right
	if (i == nx - 1)
		boundary_info.add_side(elem, 2, 2);
// Top
	if (j == ny - 1)
		boundary_info.add_side(elem, 3, 3);
// Left
	if (i == 0)
		boundary_info.add_side(elem, 4, 4);
// Front
	if (k == nz - 1)
		boundary_info.add_side(elem, 5, 5);
}

void set_boundary_names_Edge2(BoundaryInfo & boundary_info) {
	boundary_info.sideset_name(0) = "right";
	boundary_info.sideset_name(1) = "left";
}

void set_boundary_names_Quad4(BoundaryInfo & boundary_info) {
	boundary_info.sideset_name(0) = "bottom";
	boundary_info.sideset_name(1) = "right";
	boundary_info.sideset_name(2) = "top";
	boundary_info.sideset_name(3) = "left";
}

void set_boundary_names_Hex8(BoundaryInfo & boundary_info) {
	boundary_info.sideset_name(0) = "back";
	boundary_info.sideset_name(1) = "bottom";
	boundary_info.sideset_name(2) = "right";
	boundary_info.sideset_name(3) = "top";
	boundary_info.sideset_name(4) = "left";
	boundary_info.sideset_name(5) = "fron";
}

void add_node_Edge2(dof_id_type i, processor_id_type pid,
		ElemType type, MeshBase & mesh, XolotlInterface & interface) {
// Get the geometry of the Xolotl grid information
	double hy = 0.0, hz = 0.0;
	auto xolotlGrid = interface.getGridInfo(hy, hz);

// Bottom Left Back
	auto node0_ptr = mesh.add_point(libMesh::Point(xolotlGrid[i], 0.0, 0.0),
			node_id_Edge2(type, i));
	node0_ptr->set_unique_id() = node_id_Edge2(type, i);
	node0_ptr->processor_id() = pid;
}

void add_node_Quad4(dof_id_type nx, dof_id_type i, dof_id_type j,
		processor_id_type pid, ElemType type, MeshBase & mesh,
		XolotlInterface & interface) {
// Get the geometry of the Xolotl grid information
	double hy = 0.0, hz = 0.0;
	auto xolotlGrid = interface.getGridInfo(hy, hz);

// Bottom Left Back
	auto node0_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j) * hy, 0.0),
			node_id_Quad4(type, nx, i, j));
	node0_ptr->set_unique_id() = node_id_Quad4(type, nx, i, j);
	node0_ptr->processor_id() = pid;
}

void add_node_Hex8(dof_id_type nx, dof_id_type ny,
		dof_id_type i, dof_id_type j, dof_id_type k, processor_id_type pid,
		ElemType type, MeshBase & mesh, XolotlInterface & interface) {
// Get the geometry of the Xolotl grid information
	double hy = 0.0, hz = 0.0;
	auto xolotlGrid = interface.getGridInfo(hy, hz);

// Bottom Left Back
	auto node0_ptr = mesh.add_point(
			libMesh::Point(xolotlGrid[i], static_cast<Real>(j) * hy,
					static_cast<Real>(k) * hz),
			node_id_Hex8(type, nx, ny, i, j, k));
	node0_ptr->set_unique_id() = node_id_Hex8(type, nx, ny, i, j, k);
	node0_ptr->processor_id() = pid;
}

void build_cube_Edge2(UnstructuredMesh & mesh, DM da, const ElemType type,
		XolotlInterface & interface) {
	const auto pid = mesh.comm().rank();

	BoundaryInfo & boundary_info = mesh.get_boundary_info();
// xs: start grid point (not element) index on local in x direction
// xm: number of grid points owned by the local processor in x direction
// Mx: number of grid points on all processors in x direction
// xp: number of processor cores in x direction
	PetscInt xs, xm, Mx, xp;

	/* Get local grid boundaries */
	DMDAGetCorners(da, &xs, PETSC_IGNORE, PETSC_IGNORE, &xm, PETSC_IGNORE,
			PETSC_IGNORE);
	DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE, &xp,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);

	for (PetscInt i = xs; i < xs + xm; i++) {
		// We loop over grid points, but we are
		// here building elements. So that we just
		// simply skip the first x and y points since the
		// number of grid ponts is one more than
		// the number of grid elements
		if (!i)
			continue;

		dof_id_type ele_id = (i - 1);

		add_element_Edge2(da, Mx - 1, i - 1, ele_id, pid, type, mesh,
				interface);
	}

// If there is no element at the given processor
// We need to manually add all mesh nodes
	if (xs == 0 && xm == 1)
		for (PetscInt i = xs; i < xs + xm; i++)
			add_node_Edge2(i, pid, type, mesh, interface);

// Need to link up the local elements before we can know what's missing
	mesh.find_neighbors();

	mesh.find_neighbors(true);

// Set RemoteElem neighbors
	for (auto & elem_ptr : mesh.element_ptr_range())
		for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
			if (!elem_ptr->neighbor_ptr(s)
					&& !boundary_info.n_boundary_ids(elem_ptr, s))
				elem_ptr->set_neighbor(s,
						const_cast<RemoteElem *>(remote_elem));

	set_boundary_names_Edge2(boundary_info);
// Already partitioned!
	mesh.skip_partitioning(true);

// No need to renumber or find neighbors - done did it.
// Avoid deprecation message/error by _also_ setting
// allow_renumbering(false). This is a bit silly, but we want to
// catch cases where people are purely using the old "skip"
// interface and not the new flag setting one.
	mesh.allow_renumbering(false);
	mesh.prepare_for_use(/*skip_renumber (ignored!) = */false,
	/*skip_find_neighbors = */true);
}

void build_cube_Quad4(UnstructuredMesh & mesh, DM da, const ElemType type,
		XolotlInterface & interface) {
	const auto pid = mesh.comm().rank();

	BoundaryInfo & boundary_info = mesh.get_boundary_info();
// xs: start grid point (not element) index on local in x direction
// ys: start grid point index on local in y direction
// xm: number of grid points owned by the local processor in x direction
// ym: number of grid points owned by the local processor in y direction
// Mx: number of grid points on all processors in x direction
// My: number of grid points on all processors in y direction
// xp: number of processor cores in x direction
// yp: number of processor cores in y direction
	PetscInt xs, ys, xm, ym, Mx, My, xp, yp;

	/* Get local grid boundaries */
	DMDAGetCorners(da, &xs, &ys, PETSC_IGNORE, &xm, &ym, PETSC_IGNORE);
	DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE, &xp, &yp,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE);

	for (PetscInt j = ys; j < ys + ym; j++)
		for (PetscInt i = xs; i < xs + xm; i++) {
			// We loop over grid points, but we are
			// here building elements. So that we just
			// simply skip the first x and y points since the
			// number of grid ponts is one more than
			// the number of grid elements
			if (!i || !j)
				continue;

			dof_id_type ele_id = (i - 1) + (j - 1) * (Mx - 1);

			add_element_Quad4(da, Mx - 1, My - 1, i - 1, j - 1, ele_id, pid,
					type, mesh, interface);
		}

// If there is no element at the given processor
// We need to manually add all mesh nodes
	if ((ys == 0 && ym == 1) || (xs == 0 && xm == 1))
		for (PetscInt j = ys; j < ys + ym; j++)
			for (PetscInt i = xs; i < xs + xm; i++)
				add_node_Quad4(Mx, i, j, pid, type, mesh, interface);

// Need to link up the local elements before we can know what's missing
	mesh.find_neighbors();

	mesh.find_neighbors(true);

// Set RemoteElem neighbors
	for (auto & elem_ptr : mesh.element_ptr_range())
		for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
			if (!elem_ptr->neighbor_ptr(s)
					&& !boundary_info.n_boundary_ids(elem_ptr, s))
				elem_ptr->set_neighbor(s,
						const_cast<RemoteElem *>(remote_elem));

	set_boundary_names_Quad4(boundary_info);
// Already partitioned!
	mesh.skip_partitioning(true);

// No need to renumber or find neighbors - done did it.
// Avoid deprecation message/error by _also_ setting
// allow_renumbering(false). This is a bit silly, but we want to
// catch cases where people are purely using the old "skip"
// interface and not the new flag setting one.
	mesh.allow_renumbering(false);
	mesh.prepare_for_use(/*skip_renumber (ignored!) = */false,
	/*skip_find_neighbors = */true);
}

void build_cube_Hex8(UnstructuredMesh & mesh, DM da, const ElemType type,
		XolotlInterface & interface) {
	const auto pid = mesh.comm().rank();

	BoundaryInfo & boundary_info = mesh.get_boundary_info();
// xs: start grid point (not element) index on local in x direction
// ys: start grid point index on local in y direction
// zs: start grid point index on local in z direction
// xm: number of grid points owned by the local processor in x direction
// ym: number of grid points owned by the local processor in y direction
// zm: number of grid points owned by the local processor in z direction
// Mx: number of grid points on all processors in x direction
// My: number of grid points on all processors in y direction
// Mz: number of grid points on all processors in z direction
// xp: number of processor cores in x direction
// yp: number of processor cores in y direction
// zp: number of processor cores in z direction
	PetscInt xs, ys, zs, xm, ym, zm, Mx, My, Mz, xp, yp, zp;

	/* Get local grid boundaries */
	DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
	DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, &xp, &yp, &zp, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);

	for (PetscInt k = zs; k < zs + zm; k++)
		for (PetscInt j = ys; j < ys + ym; j++)
			for (PetscInt i = xs; i < xs + xm; i++) {
				// We loop over grid points, but we are
				// here building elements. So that we just
				// simply skip the first x and y points since the
				// number of grid ponts is one more than
				// the number of grid elements
				if (!i || !j || !k)
					continue;

				dof_id_type ele_id = (i - 1)
						+ (j - 1 + (k - 1) * (My - 1)) * (Mx - 1);

				add_element_Hex8(da, Mx - 1, My - 1, Mz - 1, i - 1, j - 1,
						k - 1, ele_id, pid, type, mesh, interface);
			}

// If there is no element at the given processor
// We need to manually add all mesh nodes
	if ((zs == 0 && zm == 1) || (ys == 0 && ym == 1) || (xs == 0 && xm == 1))
		for (PetscInt k = zs; k < zs + zm; k++)
			for (PetscInt j = ys; j < ys + ym; j++)
				for (PetscInt i = xs; i < xs + xm; i++)
					add_node_Hex8(Mx, My, i, j, k, pid, type, mesh, interface);

// Need to link up the local elements before we can know what's missing
	mesh.find_neighbors();

	mesh.find_neighbors(true);

// Set RemoteElem neighbors
	for (auto & elem_ptr : mesh.element_ptr_range())
		for (unsigned int s = 0; s < elem_ptr->n_sides(); s++)
			if (!elem_ptr->neighbor_ptr(s)
					&& !boundary_info.n_boundary_ids(elem_ptr, s))
				elem_ptr->set_neighbor(s,
						const_cast<RemoteElem *>(remote_elem));

	set_boundary_names_Hex8(boundary_info);
// Already partitioned!
	mesh.skip_partitioning(true);

// No need to renumber or find neighbors - done did it.
// Avoid deprecation message/error by _also_ setting
// allow_renumbering(false). This is a bit silly, but we want to
// catch cases where people are purely using the old "skip"
// interface and not the new flag setting one.
	mesh.allow_renumbering(false);
	mesh.prepare_for_use(/*skip_renumber (ignored!) = */false,
	/*skip_find_neighbors = */true);
}

void PETScDMDAMesh::buildMesh() {
// Reference to the libmesh mesh
	MeshBase & mesh = getMesh();

	MooseEnum elem_type_enum = getParam < MooseEnum > ("elem_type");

	if (!isParamValid("elem_type")) {
		// Switching on MooseEnum
		switch (_dim) {
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

	_elem_type = Utility::string_to_enum < ElemType > (elem_type_enum);

	mesh.set_mesh_dimension(_dim);
	mesh.set_spatial_dimension(_dim);

// Get the app to get the interface for the geometry of the grid
	coupling_xolotlApp * xolotl_app = dynamic_cast<coupling_xolotlApp *>(&_app);
	auto & interface = xolotl_app->getInterface();

// Switching on MooseEnum
	switch (_dim) {
	case 1:
		build_cube_Edge2(dynamic_cast<UnstructuredMesh &>(getMesh()), _dmda,
				_elem_type, interface);
		break;
	case 2:
		build_cube_Quad4(dynamic_cast<UnstructuredMesh &>(getMesh()), _dmda,
				_elem_type, interface);
		break;
        case 3:
                build_cube_Hex8(dynamic_cast<UnstructuredMesh &>(getMesh()), _dmda,
                                _elem_type, interface);
                break;
	default:
		mooseError("Does not support dimension ", _dim, "yet");
	}
}

