//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XolotlProblem.h"

registerMooseObject("coupling_xolotlApp", XolotlProblem);

template <>
InputParameters
validParams<XolotlProblem>()
{
  InputParameters params = validParams<ExternalProblem>();
  params.addRequiredParam<VariableName>("sync_variable",
                                        "The variable the solution will be synced to");
  return params;
}

XolotlProblem::XolotlProblem(const InputParameters & params)
  : ExternalProblem(params),
    _sync_to_var_name(getParam<VariableName>("sync_variable")),
    _interface(static_cast<coupling_xolotlApp &>(_app).getInterface())
{
}

void
XolotlProblem::externalSolve()
{
    // Set the time we want to reach
    _interface.setTimes(time(), dt());

    // Save the size of the dt for derivative calculation
    _dt_for_derivative = dt();

    // Save the current Xe rate
    _old_rate = *_interface.getLocalXeRate();
    
    // Run the solver
    _interface.solveXolotl();
}

void
XolotlProblem::syncSolutions(Direction direction)
{

    if (direction == Direction::FROM_EXTERNAL_APP)
    {
        auto localRate = _interface.getLocalXeRate();
        PetscInt i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
        _interface.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

        MeshBase & to_mesh = mesh().getMesh();
        auto & sync_to_var = getVariable(
                                         0, _sync_to_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);
	
        for (k = zs; k < zs + max(zm, 1); k++)
        for (j = ys; j < ys + max(ym, 1); j++)
        for (i = xs; i < xs + max(xm, 1); i++)
        {
            Node * to_node = to_mesh.node_ptr(i + (j + k * My) * Mx);
	    if (to_node->n_comp(sync_to_var.sys().number(), sync_to_var.number()) > 1)
            mooseError("Does not support multiple components");
            dof_id_type dof = to_node->dof_number(sync_to_var.sys().number(), sync_to_var.number(), 0);
	    // Compute the time derivative
	    Real value = (localRate->at(i-xs)[j-ys][k-zs] - _old_rate[i-xs][j-ys][k-zs]) / _dt_for_derivative;
//	    if (localRate->at(i-xs)[j-ys][k-zs] > 0.0) std::cout << i << " " << j << " " << localRate->at(i-xs)[j-ys][k-zs] << " " << value << std::endl;
            sync_to_var.sys().solution().set(dof, value);
        }

        sync_to_var.sys().solution().close();
    }
}
