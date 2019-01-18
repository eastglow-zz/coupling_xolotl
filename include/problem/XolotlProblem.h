//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XOLOTLPROBLEM_H
#define XOLOTLPROBLEM_H

#include "ExternalProblem.h"
#include "coupling_xolotlApp.h"
#include <interface.h>

class XolotlProblem;

template <>
InputParameters validParams<XolotlProblem>();

/**
 * This is an interface to call an external solver
  */
  class XolotlProblem : public ExternalProblem
  {
  public:
    XolotlProblem(const InputParameters & params);
    ~XolotlProblem() {}

    virtual void externalSolve() override;
    virtual void syncSolutions(Direction /*direction*/) override;

    virtual bool converged() override { return true; }

  private:
/// The name of the variable to transfer to
    const VariableName & _sync_to_var_name;
    std::string _xolotl_input_path_name;
    XolotlInterface _interface;
    std::shared_ptr<xolotlSolver::PetscSolver> _solver;
  };

#endif /* XOLOTLPROBLEM_H */
