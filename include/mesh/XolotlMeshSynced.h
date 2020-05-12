//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XOLOTLMESHSYNCED_H
#define XOLOTLMESHSYNCED_H

#include "MooseMesh.h"
#include <XolotlDLinterface.h>  // Xolotl interface

class XolotlMeshSynced;

template <>
InputParameters validParams<XolotlMeshSynced>();

/**
 * Generate a parallel (distributed) mesh from PETSc DMDA.
 * DMDA could be passed in from an application such as ExternalPetscSolverApp
 * or created on the fly. Note that this mesh object does not have one layer of
 * ghost elements. It is designed for holding the solution from an external PETSc
 * application. And then the solution can be coupled to other MOOSE-based applications
 * using the existing MultiApp transfers.
 */
class XolotlMeshSynced : public MooseMesh
{
public:
  XolotlMeshSynced(const InputParameters & parameters);
  XolotlMeshSynced(const XolotlMeshSynced & /* other_mesh */) = default;

  ~XolotlMeshSynced()
  {
  }

  // No copy
  XolotlMeshSynced & operator=(const XolotlMeshSynced & other_mesh) = delete;

  virtual void buildMesh() override;
  virtual Real getMinInDimension(unsigned int component) const override;
  virtual Real getMaxInDimension(unsigned int component) const override;
  virtual std::unique_ptr<MooseMesh> safeClone() const override;

private:
  virtual dof_id_type node_id_Edge2(const ElemType /*type*/,
                const dof_id_type /*nx*/,
                const dof_id_type /*ny*/,
                const dof_id_type i,
                const dof_id_type /*j*/,
                const dof_id_type /*k*/) const;

  virtual dof_id_type node_id_Quad4(const ElemType /*type*/,
                const dof_id_type nx,
                const dof_id_type /*ny*/,
                const dof_id_type i,
                const dof_id_type j,
                const dof_id_type /*k*/) const;

  virtual dof_id_type node_id_Hex8(const ElemType /*type*/,
                const dof_id_type nx,
                const dof_id_type ny,
                const dof_id_type nz,
                const dof_id_type i,
                const dof_id_type j,
                const dof_id_type k) const;

  virtual void add_element_Edge2(
                    const dof_id_type nx,
                    const dof_id_type i,
                    const dof_id_type elem_id,
                    const processor_id_type pid,
                    const ElemType type,
                    MeshBase & mesh) const;

  virtual void add_element_Quad4(
                    const dof_id_type nx,
                    const dof_id_type ny,
                    const dof_id_type i,
                    const dof_id_type j,
                    const dof_id_type elem_id,
                    const processor_id_type pid,
                    const ElemType type,
                    MeshBase & mesh) const;

  virtual void add_element_Hex8(
                    const dof_id_type nx,
                    const dof_id_type ny,
                    const dof_id_type nz,
                    const dof_id_type i,
                    const dof_id_type j,
                    const dof_id_type k,
                    const dof_id_type elem_id,
                    const processor_id_type pid,
                    const ElemType type,
                    MeshBase & mesh) const;

  virtual void set_boundary_names_Edge2(BoundaryInfo & boundary_info) const;
  virtual void set_boundary_names_Quad4(BoundaryInfo & boundary_info) const;
  virtual void set_boundary_names_Hex8(BoundaryInfo & boundary_info) const;
  virtual void get_indices_Edge2(const dof_id_type nx,
                    const dof_id_type /*ny*/,
                    const dof_id_type elem_id,
                    dof_id_type & i) const;

  virtual void get_indices_Quad4(const dof_id_type nx,
                    const dof_id_type /*ny*/,
                    const dof_id_type elem_id,
                    dof_id_type & i,
                    dof_id_type & j) const;

  virtual void get_indices_Hex8(const dof_id_type nx,
                    const dof_id_type ny,
                    const dof_id_type elem_id,
                    dof_id_type & i,
                    dof_id_type & j,
                    dof_id_type & k) const;

  virtual dof_id_type elem_id_Edge2(const dof_id_type i) const;

  virtual dof_id_type elem_id_Quad4(const dof_id_type nx,
                const dof_id_type /*nx*/,
                const dof_id_type i,
                const dof_id_type j,
                const dof_id_type /*k*/) const;

  virtual dof_id_type elem_id_Hex8(const dof_id_type nx,
                const dof_id_type ny,
                const dof_id_type i,
                const dof_id_type j,
                const dof_id_type k) const;

  virtual void get_neighbors_Edge2(const dof_id_type nx,
                      const dof_id_type i,
                      std::vector<dof_id_type> & neighbors) const;

  virtual void get_neighbors_Quad4(const dof_id_type nx,
                      const dof_id_type ny,
                      const dof_id_type i,
                      const dof_id_type j,
                      std::vector<dof_id_type> & neighbors) const;

  virtual void get_neighbors_Hex8(const dof_id_type nx,
                      const dof_id_type ny,
                      const dof_id_type nz,
                      const dof_id_type i,
                      const dof_id_type j,
                      const dof_id_type k,
                      std::vector<dof_id_type> & neighbors) const;

  virtual void get_ghost_neighbors_Edge2(const dof_id_type nx,
                            const MeshBase & mesh,
                            std::set<dof_id_type> & ghost_elems) const;

  virtual void get_ghost_neighbors_Quad4(const dof_id_type nx,
                      const dof_id_type ny,
                      const MeshBase & mesh,
                      std::set<dof_id_type> & ghost_elems) const;

  virtual void get_ghost_neighbors_Hex8(const dof_id_type nx,
                            const dof_id_type ny,
                            const dof_id_type nz,
                            const MeshBase & mesh,
                            std::set<dof_id_type> & ghost_elems)  const;

  virtual void add_node_Edg2(dof_id_type nx,
                dof_id_type i,
                processor_id_type pid,
                ElemType type,
                MeshBase & mesh) const;

  virtual void add_node_Qua4(dof_id_type nx,
                dof_id_type ny,
                dof_id_type i,
                dof_id_type j,
                processor_id_type pid,
                ElemType type,
                MeshBase & mesh) const;

  virtual void add_node_Hex8(dof_id_type nx,
                dof_id_type ny,
                dof_id_type nz,
                dof_id_type i,
                dof_id_type j,
                dof_id_type k,
                processor_id_type pid,
                ElemType type,
                MeshBase & mesh) const;

  virtual void build_segment_Edge2(UnstructuredMesh & mesh, const ElemType type) const;
  virtual void build_square_Quad4(UnstructuredMesh & mesh, const ElemType type) const;
  virtual void build_cube_Hex8(UnstructuredMesh & mesh, const ElemType type) const;

  virtual std::vector<double> build_xolotl_axis(int nsize, double dl) const;
  virtual int** init_xolotl_local_index_table(int ncolumn) const;
  virtual void print_xolotl_local_index_table(int **table, int ncol) const;
  virtual void fillout_xolotl_local_index_table(int** table, int ncolumn) const;
  virtual void count_num_partition_along_each_xyz(int *xp, int *yp, int *zp, int **local_index_table) const;
  virtual void get_current_process_coord(int *xpid, int *ypid, int *zpid, int pid) const;

  /// The dimension of the mesh
  MooseEnum _dim;

  /// Number of elements in x, y, z direction
  dof_id_type _nx, _ny, _nz;

  /// The min/max values for x,y,z component
  Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  /// The type of element to build
  ElemType _elem_type;

  int _moose_rank;

  /*========= Variables for the External Application =========*/
  /// Finite Difference grid parameters for Xolotl
  int _xolotl_dim;
  bool _xolotl_regulargrid;
  int _xolotl_nx, _xolotl_ny, _xolotl_nz; // Total # of grid points along each axis
  double _xolotl_dx, _xolotl_dy, _xolotl_dz; // Grid spacings
  // double *_xolotl_xc, *_xolotl_yc, *_xolotl_zc;
  std::vector<double> _xolotl_xc, _xolotl_yc, _xolotl_zc;
  std::vector<double> _xolotl_xcNR; //Non-regular grid
  double _xolotl_lx, _xolotl_ly, _xolotl_lz; // Total length of the domain along each axis
  double _xolotl_lxNR;
  int _xolotl_xi_lb, _xolotl_xi_ub; // Lower & upper bounds of the x-grid index of the MPI process
  int _xolotl_yi_lb, _xolotl_yi_ub; // Lower & upper bounds of the y-grid index of the MPI process
  int _xolotl_zi_lb, _xolotl_zi_ub; // Lower & upper bounds of the z-grid index of the MPI process
  int _xolotl_localNx, _xolotl_localNy, _xolotl_localNz;

  int _xnp, _ynp, _znp;
  int _xpid, _ypid, _zpid;

  /*
    The table of local index bounds of each Xolotl local grid at each MPI process
    Its size should be np * 6; row: mpirank & column: xs, xm, ys, ym, zs, zm
  */
  int **_xolotl_local_index_table;

  std::string _ext_lib_path_name; // External dynamic library path variable
  std::string _xolotl_input_path_name;
  void* _ext_lib_handle; // External dynamic library handle variable
  XolotlDLinterface* _xolotl_interface;
  std::shared_ptr<xolotlSolver::PetscSolver> _xolotl_solver;
};

#endif /* XOLOTLMESHSYNCED_H */
