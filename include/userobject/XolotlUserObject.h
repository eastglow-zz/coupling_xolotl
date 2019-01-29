//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XOLOTLUSEROBJECT_H
#define XOLOTLUSEROBJECT_H

#include "NodalUserObject.h"
#include "MooseVariableInterface.h"
#include <XolotlDLinterface.h>  // Xolotl interface

// Forward Declarations
class XolotlUserObject;

template <>
InputParameters validParams<NodalUserObject>();

/**
 * This UserObject computes averages of a variable storing partial
 * sums for the specified number of intervals in a direction (x,y,z).
 */
class XolotlUserObject : public NodalUserObject,
                            public MooseVariableInterface<Real>
{
public:
  XolotlUserObject(const InputParameters & parameters);

  virtual void initialize() override; //contains parameter initialization for External App.
  virtual void execute() override; //contains data transfer from MOOSE to Exernal App.
  virtual void finalize() override; //contains execution of the External App.
  virtual void threadJoin(const UserObject & y) override;
  virtual Real spatialValue(const libMesh::Point & /*p*/) const override {return calc_spatial_value_glob();} //contains data transfer from External App. to MOOSE

private:
  Real calc_spatial_value() const;
  Real calc_spatial_value_glob() const;
  virtual void map_MOOSE2Xolotl(int *rank, int *i, int *j, int *k, const Node & MOOSEnode) const;
  virtual void map_MOOSE2XolotlGlob(int *i, int *j, int *k, const Node & MOOSEnode) const;
  virtual void print_mesh_params() const;
  virtual int** init_xolotl_local_index_table(int ncolumn) const;
  virtual void print_xolotl_local_index_table(int** table, int ncolumn) const;
  virtual void fillout_xolotl_local_index_table(int** table, int ncolumn) const;
  // virtual double* build_xolotl_axis(int nsize, double dl) const;
  virtual std::vector<double> build_xolotl_axis(int nsize, double dl) const;
  virtual int max3int(int a, int b, int c) const;
  virtual int** init_xolotl_rankpair() const;
  virtual void fillout_xolotl_rankpair(int **table, int myrank, int recvRank) const;
  virtual void print_xolotl_rankpair(int **table) const;
  virtual double* vectorized_xolotl_XeRate(std::shared_ptr<xolotlSolver::PetscSolver> solver) const;
  virtual double* vectorized_xolotl_XeConc(std::shared_ptr<xolotlSolver::PetscSolver> solver) const;
  virtual double* allocate_xolotlGlobalData() const;
  virtual void localFill_xolotlGlobalXeRate(double *arr, std::shared_ptr<xolotlSolver::PetscSolver> solver) const;
  virtual void localFill_xolotlGlobalXeConc(double *arr, std::shared_ptr<xolotlSolver::PetscSolver> solver) const;
  virtual void globalFill_xolotlGlobalData(double *arr) const;
  virtual std::vector<std::tuple<int, int, int>> get_GlobalGBList(std::vector<std::tuple<int, int, int>> gbLocal) const;
  virtual int iiGlob(int i, int j, int k) const;
  virtual int ii(int i, int j, int k) const;
  virtual int iiR(int rank, int i, int j, int k) const;

  MooseVariable & _var;
  const VariableValue & _u;
  const VariableValue & _v;
  int _moose_rank;
  std::vector<std::tuple<int, int, int>> _GBListLocal;
  std::vector<std::tuple<int, int, int>> _GBList;

  /*========= Variables for the External Application =========*/
  /// Finite Difference grid parameters for Xolotl
  int _xolotl_dim;
  bool _xolotl_regulargrid;
  int _xolotl_nx, _xolotl_ny, _xolotl_nz; // Total # of grid points along each axis
  double _xolotl_dx, _xolotl_dy, _xolotl_dz; // Grid spacings
  // double *_xolotl_xc, *_xolotl_yc, *_xolotl_zc;
  std::vector<double> _xolotl_xc, _xolotl_yc, _xolotl_zc;
  double _xolotl_lx, _xolotl_ly, _xolotl_lz; // Total length of the domain along each axis
  int _xolotl_xi_lb, _xolotl_xi_ub; // Lower & upper bounds of the x-grid index of the MPI process
  int _xolotl_yi_lb, _xolotl_yi_ub; // Lower & upper bounds of the y-grid index of the MPI process
  int _xolotl_zi_lb, _xolotl_zi_ub; // Lower & upper bounds of the z-grid index of the MPI process
  int _xolotl_localNx, _xolotl_localNy, _xolotl_localNz;

  /*
    The table of local index bounds of each Xolotl local grid at each MPI process
    Its size should be np * 6; row: mpirank & column: xs, xm, ys, ym, zs, zm
  */
  int **_xolotl_local_index_table;

  /*
    The table of pairs of ranks for send to or recv from
    Its size should be np; row: _moose_rank & coulmn: sendRank, recvRank
  */
  int **_xolotl_rank_pair;

  std::string _ext_lib_path_name; // External dynamic library path variable
  std::string _xolotl_input_path_name;
  void* _ext_lib_handle; // External dynamic library handle variable
  XolotlDLinterface* _xolotl_interface;
  int _argc = 3;
  char ** _argv = new char*[_argc];
  std::string _parameterFile = "crap";
  // std::shared_ptr<xolotlCore::Options> _xolotl_options;
  std::shared_ptr<xolotlSolver::PetscSolver> _xolotl_solver;
  double *_xolotl_XeRate;
  double *_xolotl_XeConc;

  double *_xolotl_GlobalXeRate;
  double *_xolotl_GlobalXeConc;
};

#endif //MYDIFFUSIONUSEROBJECT_H
