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

template<>
InputParameters validParams<XolotlProblem>() {
	InputParameters params = validParams<ExternalProblem>();
	params.addRequiredParam < VariableName
			> ("sync_variable", "The variable the solution will be synced to");
	params.addRequiredParam < VariableName
			> ("sync_GB", "The variable the GB will be synced to");
	return params;
}

XolotlProblem::XolotlProblem(const InputParameters & params) :
		ExternalProblem(params), _sync_to_var_name(
				getParam < VariableName > ("sync_variable")), _sync_from_var_name(
				getParam < VariableName > ("sync_GB")), _interface(
				static_cast<coupling_xolotlApp &>(_app).getInterface()) {
	PetscInt xs, ys, zs, xm, ym, zm, Mx, My, Mz;
        _interface.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
	
	// Initialize old rate
	_old_rate.clear();
	for (int i = 0; i < max(xm, 1); i++) {
		std::vector<std::vector<double> > tempTempVector;
		for (int j = 0; j < max(ym, 1); j++) {
			std::vector<double> tempVector;
			for (int k = 0; k < max(zm, 1); k++) {
				tempVector.push_back(0.0);
			}
			tempTempVector.push_back(tempVector);
		}
		_old_rate.push_back(tempTempVector);
	}
}

void XolotlProblem::externalSolve() {
	// Make sure the times will be the same for every MPI rank
//	Real xolotlTime = 0.0, xolotlDT = 0.0;
//	MPI_Allreduce(&time(), &xolotlTime, 1, MPI_DOUBLE, MPI_MIN, _communicator.get());
//	MPI_Allreduce(&dt(), &xolotlDT, 1, MPI_DOUBLE, MPI_MIN, _communicator.get());

	// Set the time we want to reach
	_interface.setTimes(time(), dt());
	// Reset the concentrations where the GBs are
	_interface.initGBLocation();
	// Save the size of the dt for derivative calculation
	_dt_for_derivative = dt();
	// Run the solver
	_interface.solveXolotl();
}

void XolotlProblem::syncSolutions(Direction direction) {
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_interface.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
	if (direction == Direction::FROM_EXTERNAL_APP) {
		MeshBase & to_mesh = mesh().getMesh();
		auto & sync_to_var = getVariable(0, _sync_to_var_name,
				Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);

		for (k = zs; k < zs + max(zm, 1); k++)
			for (j = ys; j < ys + max(ym, 1); j++)
				for (i = xs; i < xs + max(xm, 1); i++) {
					Node * to_node = to_mesh.node_ptr(i + (j + k * My) * Mx);
					if (to_node->n_comp(sync_to_var.sys().number(),
							sync_to_var.number()) > 1)
						mooseError("Does not support multiple components");
					dof_id_type dof = to_node->dof_number(
							sync_to_var.sys().number(), sync_to_var.number(),
							0);
					// Compute the time derivative
					Real current_rate = _interface.getLocalXeRate(i - xs, j - ys, k - zs);
					Real value = (current_rate
							- _old_rate[i - xs][j - ys][k - zs])
							/ _dt_for_derivative;
					sync_to_var.sys().solution().set(dof, value);
					// Update the old rate
					_old_rate[i - xs][j - ys][k - zs] = current_rate;
				}

		sync_to_var.sys().solution().close();
	}

	if (direction == Direction::TO_EXTERNAL_APP) {
		MeshBase & to_mesh = mesh().getMesh();
		auto & sync_from_var = getVariable(0, _sync_from_var_name,
				Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);

		// Create a list of GB
		std::vector<int> localGBList;

		for (k = zs; k < zs + max(zm, 1); k++)
			for (j = ys; j < ys + max(ym, 1); j++)
				for (i = xs; i < xs + max(xm, 1); i++) {
					Node * to_node = to_mesh.node_ptr(i + (j + k * My) * Mx);
					if (to_node->n_comp(sync_from_var.sys().number(),
							sync_from_var.number()) > 1)
						mooseError("Does not support multiple components");
					dof_id_type dof = to_node->dof_number(
							sync_from_var.sys().number(),
							sync_from_var.number(), 0);
					// Get the value
					Real value = sync_from_var.sys().solution()(dof);
					// Test if it is a GB
					if (value < 0.9) {
						localGBList.push_back(i);
						localGBList.push_back(j);
						localGBList.push_back(k);
					}
				}

		sync_from_var.sys().solution().close();

		// Clear the GB list
		_gb_list.clear();
		// Prepare for all gather v
		int worldSize = 0;
		MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
		std::vector<int> sizes;
		for (int m = 0; m < worldSize; m++) {
			sizes.push_back(0);
		}
		int n = localGBList.size(); 
		MPI_Allgather(&n, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
		std::vector<int> offsets;
		int s = 0;
		for (int m = 0; m < worldSize; m++) {
			offsets.push_back(s);
			s += sizes[m];
		}
		int N = s;
		_gb_list.resize(N);
		_interface.resetGBVector();
		// Gather all the GB locations on all the procs
		MPI_Allgatherv((void*) localGBList.data(), localGBList.size(), MPI_INT,
				(void*) _gb_list.data(), sizes.data(), offsets.data(), MPI_INT,
				MPI_COMM_WORLD);
		for (int m = 0; m < _gb_list.size(); m += 3) {
			_interface.setGBLocation(_gb_list[m], _gb_list[m + 1], _gb_list[m + 2]);
		}
	}
}
