#include "coupling_xolotlApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "XolotlProblem.h"
#include "Executioner.h"
#include "ModulesApp.h"

InputParameters coupling_xolotlApp::validParams() {
	InputParameters params = MooseApp::validParams();

	return params;
}

coupling_xolotlApp::coupling_xolotlApp(InputParameters parameters) :
		MooseApp(parameters) {
	coupling_xolotlApp::registerAll(_factory, _action_factory, _syntax);
}

coupling_xolotlApp::~coupling_xolotlApp() {
}

void coupling_xolotlApp::registerAll(Factory & f, ActionFactory & af,
		Syntax & s) {
	ModulesApp::registerAll(f, af, s);
	Registry::registerObjectsTo(f, { "coupling_xolotlApp" });
	Registry::registerActionsTo(af, { "coupling_xolotlApp" });

	/* register custom execute flags, action syntax, etc. here */
}

void coupling_xolotlApp::registerApps() {
	registerApp (coupling_xolotlApp);
}

std::shared_ptr<Backup> coupling_xolotlApp::backup() {
	// Get the state from Xolotl
	mooseAssert(_executioner, "Executioner is nullptr");
	XolotlProblem & xolotl_problem = (XolotlProblem &) _executioner->feProblem();
	xolotl_problem.saveState();

	// Back it up
	return MooseApp::backup();
}

void coupling_xolotlApp::restore(std::shared_ptr<Backup> backup,
		bool for_restart) {
	// Restore the state
	MooseApp::restore(backup, for_restart);

	// Set it in Xolotl
	mooseAssert(_executioner, "Executioner is nullptr");
        XolotlProblem & xolotl_problem = (XolotlProblem &) _executioner->feProblem();
	xolotl_problem.setState();
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void coupling_xolotlApp__registerAll(Factory & f, ActionFactory & af,
		Syntax & s) {
	coupling_xolotlApp::registerAll(f, af, s);
}

extern "C" void coupling_xolotlApp__registerApps() {
	coupling_xolotlApp::registerApps();
}
