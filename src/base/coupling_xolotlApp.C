#include "coupling_xolotlApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<coupling_xolotlApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

coupling_xolotlApp::coupling_xolotlApp(InputParameters parameters) : MooseApp(parameters)
{
  coupling_xolotlApp::registerAll(_factory, _action_factory, _syntax);
}

coupling_xolotlApp::~coupling_xolotlApp() {}

void
coupling_xolotlApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"coupling_xolotlApp"});
  Registry::registerActionsTo(af, {"coupling_xolotlApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
coupling_xolotlApp::registerApps()
{
  registerApp(coupling_xolotlApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
coupling_xolotlApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  coupling_xolotlApp::registerAll(f, af, s);
}

extern "C" void
coupling_xolotlApp__registerApps()
{
  coupling_xolotlApp::registerApps();
}
