//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "coupling_xolotlTestApp.h"
#include "coupling_xolotlApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<coupling_xolotlTestApp>()
{
  InputParameters params = validParams<coupling_xolotlApp>();
  return params;
}

coupling_xolotlTestApp::coupling_xolotlTestApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  coupling_xolotlApp::registerObjectDepends(_factory);
  coupling_xolotlApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  coupling_xolotlApp::associateSyntaxDepends(_syntax, _action_factory);
  coupling_xolotlApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  coupling_xolotlApp::registerExecFlags(_factory);

  bool use_test_objs = getParam<bool>("allow_test_objects");
  if (use_test_objs)
  {
    coupling_xolotlTestApp::registerObjects(_factory);
    coupling_xolotlTestApp::associateSyntax(_syntax, _action_factory);
    coupling_xolotlTestApp::registerExecFlags(_factory);
  }
}

coupling_xolotlTestApp::~coupling_xolotlTestApp() {}

void
coupling_xolotlTestApp::registerApps()
{
  registerApp(coupling_xolotlApp);
  registerApp(coupling_xolotlTestApp);
}

void
coupling_xolotlTestApp::registerObjects(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new test objects here! */
}

void
coupling_xolotlTestApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
  /* Uncomment Syntax and ActionFactory parameters and register your new test objects here! */
}

void
coupling_xolotlTestApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execute flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
coupling_xolotlTestApp__registerApps()
{
  coupling_xolotlTestApp::registerApps();
}

// External entry point for dynamic object registration
extern "C" void
coupling_xolotlTestApp__registerObjects(Factory & factory)
{
  coupling_xolotlTestApp::registerObjects(factory);
}

// External entry point for dynamic syntax association
extern "C" void
coupling_xolotlTestApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  coupling_xolotlTestApp::associateSyntax(syntax, action_factory);
}

// External entry point for dynamic execute flag loading
extern "C" void
coupling_xolotlTestApp__registerExecFlags(Factory & factory)
{
  coupling_xolotlTestApp::registerExecFlags(factory);
}
