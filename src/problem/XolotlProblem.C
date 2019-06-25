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
				static_cast<coupling_xolotlApp &>(_app).getInterface()), _old_rate(
				declareRestartableData
						< std::vector<std::vector<std::vector<Real> > >
						> ("old_rate")), _current_time(
				declareRestartableData < Real > ("current_time", 0.0)), _current_dt(
				declareRestartableData < Real > ("current_dt", 0.0)), _previous_time(
				declareRestartableData < Real > ("previous_time", 0.0)), _n_xenon(
				declareRestartableData < Real > ("n_xenon", 0.0)), _previous_xe_flux(
				declareRestartableData
						< std::vector<std::vector<std::vector<Real> > >
						> ("previous_xe_flux")), _local_rate(
				declareRestartableData
						< std::vector<std::vector<std::vector<Real> > >
						> ("local_rate")), _conc_vector(
				declareRestartableData
						< std::vector<
								std::vector<
										std::vector<
												std::vector<std::pair<int, Real> > > > >
						> ("conc_vector")) {
	PetscInt xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_interface.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Initialize old rate
	_old_rate.clear();
	for (int i = 0; i < max(xm, 1); i++) {
		std::vector < std::vector<double> > tempTempVector;
		for (int j = 0; j < max(ym, 1); j++) {
			std::vector<double> tempVector;
			for (int k = 0; k < max(zm, 1); k++) {
				tempVector.push_back(0.0);
			}
			tempTempVector.push_back(tempVector);
		}
		_old_rate.push_back(tempTempVector);
	}

	// Initialize has run for Xolotl
	_xolotl_has_run = false;
}

void XolotlProblem::externalSolve() {
	_xolotl_has_run = false;
	// Check that the next time is larger than the current one
	if (time() > _current_time) {
		// Set the time we want to reach
		_interface.setTimes(time(), dt());
		// Reset the concentrations where the GBs are
		_interface.initGBLocation();
		// Save the size of the dt for derivative calculation
		_dt_for_derivative = dt();
		// Run the solver
		_interface.solveXolotl();
		// Save the current time
		_current_time = time();
		// Set Xolotl has run
		_xolotl_has_run = true;
	}
}

void XolotlProblem::syncSolutions(Direction direction) {
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_interface.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
	if (direction == Direction::FROM_EXTERNAL_APP && _xolotl_has_run) {
		MeshBase & to_mesh = mesh().getMesh();
		auto & sync_to_var = getVariable(0, _sync_to_var_name,
				Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);

		// Get the rate vector
		auto rate_vector = _interface.getLocalXeRate();

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
					Real current_rate = rate_vector[i - xs][j - ys][k - zs];
					Real value = (current_rate
							- _old_rate[i - xs][j - ys][k - zs])
							/ _dt_for_derivative;
					sync_to_var.sys().solution().set(dof, value);
				}

		// Update the old rate
		_old_rate = rate_vector;

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

		// Clear the GB list in Xolotl
		_interface.resetGBVector();
		// Pass the local list to Xolotl
		for (int m = 0; m < localGBList.size(); m += 3) {
			_interface.setGBLocation(localGBList[m], localGBList[m + 1],
					localGBList[m + 2]);
		}
	}
}

void XolotlProblem::saveState() {
	// Update the values from Xolotl
	_conc_vector = _interface.getConcVector();
	_local_rate = _interface.getLocalXeRate();
	_previous_xe_flux = _interface.getPreviousXeFlux();
	_current_dt = _interface.getCurrentDt();
	_previous_time = _interface.getPreviousTime();
	_n_xenon = _interface.getNXeGB();
	_old_rate = _interface.getLocalXeRate();
}

void XolotlProblem::setState() {
	// Set them in Xolotl
	_interface.setConcVector(_conc_vector);
	_interface.setLocalXeRate(_local_rate);
	_interface.setPreviousXeFlux(_previous_xe_flux);
	_interface.setCurrentTimes(_current_time, _current_dt);
	_interface.setPreviousTime(_previous_time);
	_interface.setNXeGB(_n_xenon);
}
