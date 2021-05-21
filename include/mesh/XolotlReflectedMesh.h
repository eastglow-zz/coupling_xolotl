//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XOLOTLREFLECTEDAMESH_H
#define XOLOTLREFLECTEDAMESH_H

#include "MooseMesh.h"

class XolotlReflectedMesh;

template<>
InputParameters validParams<XolotlReflectedMesh>();

/**
 * Generate a parallel (distributed) mesh from PETSc DMDA.
 * DMDA could be passed in from an application such as ExternalPetscSolverApp
 * or created on the fly. Note that this mesh object does not have one layer of
 * ghost elements. It is designed for holding the solution from an external PETSc
 * application. And then the solution can be coupled to other MOOSE-based applications
 * using the existing MultiApp transfers.
 */
class XolotlReflectedMesh: public MooseMesh {
public:
	XolotlReflectedMesh(const InputParameters &parameters);
	XolotlReflectedMesh(const XolotlReflectedMesh& /* other_mesh */) = default;

	~XolotlReflectedMesh() {
	}

	// No copy
	XolotlReflectedMesh& operator=(const XolotlReflectedMesh &other_mesh) = delete;

	virtual std::unique_ptr<MooseMesh> safeClone() const override;

	virtual void buildMesh() override;

protected:
	/// The path to the input file for Xolotl
	FileName _xolotl_input_path_name;

	/// The dimension of the mesh
	MooseEnum _dim;

	/// Number of elements in x, y, z direction
	dof_id_type _nx, _ny, _nz;

	/// The min/max values for x,y,z component
	Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

	/// Mesh object
	DM _dmda;
};

#endif /* XOLOTLREFLECTEDMESH_H */
