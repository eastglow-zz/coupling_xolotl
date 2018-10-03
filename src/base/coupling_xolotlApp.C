#include "coupling_xolotlApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"


// My custom Executioner
#include "XolotlExecutioner.h"

// My custom TimeStepper
#include "XolotlTimeStepper.h"

template <>
InputParameters
validParams<coupling_xolotlApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

coupling_xolotlApp::coupling_xolotlApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  coupling_xolotlApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  coupling_xolotlApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  coupling_xolotlApp::registerExecFlags(_factory);
}

coupling_xolotlApp::~coupling_xolotlApp() {}

void
coupling_xolotlApp::registerApps()
{
  registerApp(coupling_xolotlApp);
}

void
coupling_xolotlApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"coupling_xolotlApp"});
    //Registering XolotlExecutioner
    registerExecutioner(XolotlExecutioner);

    //Registering XolotlTimeStepper
    registerTimeStepper(XolotlTimeStepper);
}

void
coupling_xolotlApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"coupling_xolotlApp"});

  /* Uncomment Syntax parameter and register your new production objects here! */
}

void
coupling_xolotlApp::registerObjectDepends(Factory & /*factory*/)
{
}

void
coupling_xolotlApp::associateSyntaxDepends(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

void
coupling_xolotlApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execution flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
coupling_xolotlApp__registerApps()
{
  coupling_xolotlApp::registerApps();
}

extern "C" void
coupling_xolotlApp__registerObjects(Factory & factory)
{
  coupling_xolotlApp::registerObjects(factory);
}

extern "C" void
coupling_xolotlApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  coupling_xolotlApp::associateSyntax(syntax, action_factory);
}

extern "C" void
coupling_xolotlApp__registerExecFlags(Factory & factory)
{
  coupling_xolotlApp::registerExecFlags(factory);
}
