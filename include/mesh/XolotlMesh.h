//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XOLOTLMESH_H
#define XOLOTLMESH_H

#include "MooseMesh.h"
#include <XolotlDLinterface.h>  // Xolotl interface
// #include <vector>

class XolotlMesh;

template <>
InputParameters validParams<XolotlMesh>();

/**
 * Mesh generated from parameters
 */
class XolotlMesh : public MooseMesh
{
public:
  XolotlMesh(const InputParameters & parameters);
  XolotlMesh(const XolotlMesh & /* other_mesh */) = default;

  // No copy
  XolotlMesh & operator=(const XolotlMesh & other_mesh) = delete;

  virtual std::unique_ptr<MooseMesh> safeClone() const override;

  virtual void buildMesh() override;
  virtual Real getMinInDimension(unsigned int component) const override;
  virtual Real getMaxInDimension(unsigned int component) const override;
private:
  // virtual double* build_xolotl_axis(int nsize, double dl) const;
  virtual std::vector<double> build_xolotl_axis(int nsize, double dl) const;
  virtual void map_MOOSE2XolotlGlob(int *ireturn, int *jreturn, int *kreturn, const Node & MOOSEnode) const;
  virtual int max3int(int a, int b, int c) const;

protected:
  /// The dimension of the mesh
  MooseEnum _dim;

  /// Number of elements in x, y, z direction
  unsigned int _nx, _ny, _nz;

  /// The min/max values for x,y,z component
  Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  /// All of the libmesh build_line/square/cube routines support an
  /// option to grade the mesh into the boundaries according to the
  /// spacing of the Gauss-Lobatto quadrature points.  Defaults to
  /// false, and cannot be used in conjunction with x, y, and z
  /// biasing.
  bool _gauss_lobatto_grid;

  /// The amount by which to bias the cells in the x,y,z directions.
  /// Must be in the range 0.5 <= _bias_x <= 2.0.
  /// _bias_x < 1 implies cells are shrinking in the x-direction.
  /// _bias_x==1 implies no bias (original mesh unchanged).
  /// _bias_x > 1 implies cells are growing in the x-direction.
  Real _bias_x, _bias_y, _bias_z;



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

  std::string _ext_lib_path_name; // External dynamic library path variable
  std::string _xolotl_input_path_name;
  void* _ext_lib_handle; // External dynamic library handle variable
  XolotlDLinterface* _xolotl_interface;
  std::shared_ptr<xolotlSolver::PetscSolver> _xolotl_solver;
};

#endif /* XOLOTLMESH_H */
