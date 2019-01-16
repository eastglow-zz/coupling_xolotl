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
  virtual Real spatialValue(const libMesh::Point & /*p*/) const override {return calc_spatial_value();} //contains data transfer from External App. to MOOSE

private:
  Real calc_spatial_value() const;
  virtual Real getMinInDimension(unsigned int component) const;
  virtual Real getMaxInDimension(unsigned int component) const;
  virtual std::vector<libMesh::Point> initExtCoords() const;
  virtual unsigned int map_MOOSE2Ext(const Node & MOOSEnode) const;
  virtual void print_mesh_params() const;
  virtual void print_ext_coord(std::vector<libMesh::Point> a) const;

  std::vector<libMesh::Point> _ext_coord;
  std::vector<Real> _ext_data;

  /// The dimension of the mesh
  MooseEnum _dim;

  /// Number of elements in x, y, z direction
  unsigned int _nx, _ny, _nz;

  /// The min/max values for x,y,z component (from input file)
  Real _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  MooseVariable & _var;
  const VariableValue & _u;
  const VariableValue & _v;

  /*========= Variables for the External Application =========*/
  /// Finite Difference mesh parameters for the External App.
  Real xmin, xmax, ymin, ymax, zmin, zmax;
  unsigned int nNode_x, nNode_y, nNode_z;
  Real dx, dy, dz;
  std::string _ext_lib_path_name; // External dynamic library path variable
  std::string _xolotl_input_path_name;
  void* _ext_lib_handle; // External dynamic library handle variable
  XolotlDLinterface* _xolotl_interface;
  int _argc = 3;
  char ** _argv = new char*[_argc];
  std::string _parameterFile = "crap";
  std::shared_ptr<xolotlSolver::PetscSolver> _xolotl_solver;
  std::vector<std::vector<std::vector<double> > > * _xolotl_LocalXeRate;
};

#endif //MYDIFFUSIONUSEROBJECT_H
