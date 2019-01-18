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
  params.addRequiredParam<std::string>("XolotlInput_path_name",
                                       "Name with the path for the Xolotl input file");
  return params;
}

XolotlProblem::XolotlProblem(const InputParameters & params)
  : ExternalProblem(params),
    _sync_to_var_name(getParam<VariableName>("sync_variable")),
    _xolotl_input_path_name(getParam<std::string>("XolotlInput_path_name"))
{
    int argc = 3;
    char ** argv = new char*[argc];
    std::string parameterFile = "bla";
    argv[0] = new char[parameterFile.length() + 1];
    strcpy(argv[0], parameterFile.c_str());
    argv[1] = new char[_xolotl_input_path_name.length() + 1];
    strcpy(argv[1], _xolotl_input_path_name.c_str());
    argv[2] = 0; // null-terminate the array
    
    _solver = _interface.initializeXolotl(argc,
              argv, MPI_COMM_WORLD, false);
}

void
XolotlProblem::externalSolve()
{
    // Set the time we want to reach
    _interface.setTimes(_solver, time(), dt());
    
    // Run the solver
    _interface.solveXolotl(_solver);
}

void
XolotlProblem::syncSolutions(Direction direction)
{
    if (direction == Direction::FROM_EXTERNAL_APP)
    {
        auto localRate = _interface.getLocalXeRate(_solver);
        PetscInt i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
        _interface.getLocalCoordinates(_solver, xs, xm, Mx, ys, ym, My, zs, zm, Mz);

        MeshBase & to_mesh = mesh().getMesh();
        auto & sync_to_var = getVariable(
                                         0, _sync_to_var_name, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD);

        for (k = zs; k < zs + max(zm, 1); k++)
        for (j = ys; j < ys + max(ym, 1); j++)
        for (i = xs; i < xs + max(xm, 1); i++)
        {
            Node * to_node = to_mesh.node_ptr(i + j * Mx);
            if (to_node->n_comp(sync_to_var.sys().number(), sync_to_var.number()) > 1)
            mooseError("Does not support multiple components");
            dof_id_type dof = to_node->dof_number(sync_to_var.sys().number(), sync_to_var.number(), 0);
            sync_to_var.sys().solution().set(dof, localRate->at(i-xs)[j-ys][k-zs]);
        }

        sync_to_var.sys().solution().close();
    }
}
