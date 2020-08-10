//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "coupling_xolotlTestApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters coupling_xolotlTestApp::validParams() {
	InputParameters params = coupling_xolotlApp::validParams();
	return params;
}

coupling_xolotlTestApp::coupling_xolotlTestApp(InputParameters parameters) :
		coupling_xolotlApp(parameters) {
	coupling_xolotlTestApp::registerAll(_factory, _action_factory, _syntax,
			getParam<bool>("allow_test_objects"));
}

coupling_xolotlTestApp::~coupling_xolotlTestApp() {
}

void coupling_xolotlTestApp::registerAll(Factory &f, ActionFactory &af,
		Syntax &s, bool use_test_objs) {
	coupling_xolotlApp::registerAll(f, af, s);
	if (use_test_objs) {
		Registry::registerObjectsTo(f, { "coupling_xolotlTestApp" });
		Registry::registerActionsTo(af, { "coupling_xolotlTestApp" });
	}
}

void coupling_xolotlTestApp::registerApps() {
	registerApp (coupling_xolotlApp);
	registerApp (coupling_xolotlTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void coupling_xolotlTestApp__registerAll(Factory &f,
		ActionFactory &af, Syntax &s) {
	coupling_xolotlTestApp::registerAll(f, af, s);
}
extern "C" void coupling_xolotlTestApp__registerApps() {
	coupling_xolotlTestApp::registerApps();
}
