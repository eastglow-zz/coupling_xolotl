//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include <cmath>
#include "XolotlProblem.h"
#include "SystemBase.h"

using std::max;

registerMooseObject("coupling_xolotlApp", XolotlProblem);

template<>
InputParameters validParams<XolotlProblem>() {
	InputParameters params = validParams<ExternalProblem>();
	params.addRequiredParam < VariableName
			> ("sync_rate", "The variable the rate will be synced to");
	params.addRequiredParam < VariableName
			> ("sync_GB", "The variable the GB will be synced to");
	params.addRequiredParam < VariableName
			> ("sync_mono", "The variable the monomer will be synced to");
	params.addRequiredParam < VariableName
			> ("sync_frac", "The variable the V fraction will be synced to");
	return params;
}
template<>
void dataStore(std::ostream &stream, std::tuple<Real, Real, Real, Real> &foo,
		void *context) {
	storeHelper(stream, std::get < 0 > (foo), context);
	storeHelper(stream, std::get < 1 > (foo), context);
	storeHelper(stream, std::get < 2 > (foo), context);
	storeHelper(stream, std::get < 3 > (foo), context);
}

template<>
void dataLoad(std::istream &stream, std::tuple<Real, Real, Real, Real> &foo,
		void *context) {
	loadHelper(stream, std::get < 0 > (foo), context);
	loadHelper(stream, std::get < 1 > (foo), context);
	loadHelper(stream, std::get < 2 > (foo), context);
	loadHelper(stream, std::get < 3 > (foo), context);
}

XolotlProblem::XolotlProblem(const InputParameters &params) :
		ExternalProblem(params), _sync_rate(
				getParam < VariableName > ("sync_rate")), _sync_gb(
				getParam < VariableName > ("sync_GB")), _sync_mono(
				getParam < VariableName > ("sync_mono")), _sync_frac(
				getParam < VariableName > ("sync_frac")), _interface(
				static_cast<coupling_xolotlApp&>(_app).getInterface()), _old_rate(
				declareRestartableData
						< std::vector<std::vector<std::vector<Real> > >
						> ("old_rate")), _current_time(
				declareRestartableData < Real > ("current_time", 0.0)), _current_dt(
				declareRestartableData < Real > ("current_dt", 0.0)), _previous_time(
				declareRestartableData < Real > ("previous_time", 0.0)), _n_xenon(
				declareRestartableData < Real > ("n_xenon", 0.0)), _local_NE(
				declareRestartableData
						< std::vector<
								std::vector<std::vector<std::array<Real, 4> > > >
						> ("local_NE")), _conc_vector(
				declareRestartableData
						< std::vector<
								std::vector<
										std::vector<
												std::vector<std::pair<xolotl::IdType, Real> > > > >
						> ("conc_vector")) {
	xolotl::IdType xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_interface->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);

	// Initialize old rate
	_old_rate.clear();
	for (int i = 0; i < max(xm, (xolotl::IdType) 1); i++) {
		std::vector < std::vector<double> > tempTempVector;
		for (int j = 0; j < max(ym, (xolotl::IdType) 1); j++) {
			std::vector<double> tempVector;
			for (int k = 0; k < max(zm, (xolotl::IdType) 1); k++) {
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
		_interface->setTimes(time(), dt());
		// Save the size of the dt for derivative calculation
		_dt_for_derivative = dt();
		// Run the solver
		_interface->solveXolotl();
		// Save the current time
		_current_time = time();
		// Set Xolotl has run
		_xolotl_has_run = true;
	}
}

bool XolotlProblem::converged() {
	bool conv = _interface->getConvergenceStatus();
	return conv;
}

void XolotlProblem::syncSolutions(Direction direction) {
	xolotl::IdType i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_interface->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
	if (direction == Direction::FROM_EXTERNAL_APP && _xolotl_has_run) {
		MeshBase &to_mesh = mesh().getMesh();
		auto &sync_rate = getVariable(0, _sync_rate,
				Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);
		auto &sync_mono = getVariable(0, _sync_mono,
				Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);
		auto &sync_frac = getVariable(0, _sync_frac,
				Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);

		// Get the rate vector
		auto ne_vector = _interface->getLocalNE();

		for (k = zs; k < zs + max(zm, (xolotl::IdType) 1); k++)
			for (j = ys; j < ys + max(ym, (xolotl::IdType) 1); j++)
				for (i = xs; i < xs + max(xm, (xolotl::IdType) 1); i++) {
					Node *to_node = to_mesh.node_ptr(i + (j + k * My) * Mx);
					if (to_node->n_comp(sync_rate.sys().number(),
							sync_rate.number()) > 1)
						mooseError("Does not support multiple components");
					dof_id_type dof_rate = to_node->dof_number(
							sync_rate.sys().number(), sync_rate.number(), 0);
					dof_id_type dof_mono = to_node->dof_number(
							sync_mono.sys().number(), sync_mono.number(), 0);
					dof_id_type dof_frac = to_node->dof_number(
							sync_frac.sys().number(), sync_frac.number(), 0);
					// Compute the time derivative
					Real current_rate = std::get < 0
							> (ne_vector[i - xs][j - ys][k - zs]);
					Real value = (current_rate
							- _old_rate[i - xs][j - ys][k - zs])
							/ _dt_for_derivative;
					sync_rate.sys().solution().set(dof_rate, value);
					sync_mono.sys().solution().set(dof_mono,
							std::get < 2 > (ne_vector[i - xs][j - ys][k - zs]));
					sync_frac.sys().solution().set(dof_frac,
							std::get < 3 > (ne_vector[i - xs][j - ys][k - zs]));

					// Update the old rate
					_old_rate[i - xs][j - ys][k - zs] = current_rate;
				}

		sync_rate.sys().solution().close();
		sync_mono.sys().solution().close();
		sync_frac.sys().solution().close();
	}

	if (direction == Direction::TO_EXTERNAL_APP) {
		MeshBase &to_mesh = mesh().getMesh();
		auto &sync_gb = getVariable(0, _sync_gb, Moose::VarKindType::VAR_ANY,
				Moose::VarFieldType::VAR_FIELD_STANDARD);

		// Create a list of GB
		std::vector<int> localGBList;

		for (k = zs; k < zs + max(zm, (xolotl::IdType) 1); k++)
			for (j = ys; j < ys + max(ym, (xolotl::IdType) 1); j++)
				for (i = xs; i < xs + max(xm, (xolotl::IdType) 1); i++) {
					Node *to_node = to_mesh.node_ptr(i + (j + k * My) * Mx);
					if (to_node->n_comp(sync_gb.sys().number(),
							sync_gb.number()) > 1)
						mooseError("Does not support multiple components");
					dof_id_type dof = to_node->dof_number(
							sync_gb.sys().number(), sync_gb.number(), 0);
					// Get the value
					Real value = sync_gb.sys().solution()(dof);
					// Test if it is a GB
					if (value < 0.9) {
						localGBList.push_back(i);
						localGBList.push_back(j);
						localGBList.push_back(k);
					}
				}

		sync_gb.sys().solution().close();

		// Clear the GB list in Xolotl
		_interface->resetGBVector();
		// Pass the local list to Xolotl
		for (int m = 0; m < localGBList.size(); m += 3) {
			_interface->setGBLocation(localGBList[m], localGBList[m + 1],
					localGBList[m + 2]);
		}
	}
}

void XolotlProblem::saveState() {
	// Update the values from Xolotl
	_conc_vector = _interface->getConcVector();
	_local_NE = _interface->getLocalNE();
	_current_dt = _interface->getCurrentDt();
	_previous_time = _interface->getPreviousTime();
	_n_xenon = _interface->getNXeGB();

	xolotl::IdType i, j, k, xs, ys, zs, xm, ym, zm, Mx, My, Mz;
	_interface->getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
	// Set old rate from local NE
	for (k = zs; k < zs + max(zm, (xolotl::IdType) 1); k++)
		for (j = ys; j < ys + max(ym, (xolotl::IdType) 1); j++)
			for (i = xs; i < xs + max(xm, (xolotl::IdType) 1); i++) {
				_old_rate[i - xs][j - ys][k - zs] = std::get < 0
						> (_local_NE[i - xs][j - ys][k - zs]);
			}
}

void XolotlProblem::setState() {
	// Set them in Xolotl
	_interface->setConcVector(_conc_vector);
	_interface->setLocalNE(_local_NE);
	_interface->setCurrentTimes(_current_time, _current_dt);
	_interface->setPreviousTime(_previous_time);
	_interface->setNXeGB(_n_xenon);
}
