//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef COUPLING_XOLOTLTESTAPP_H
#define COUPLING_XOLOTLTESTAPP_H

#include "MooseApp.h"

class coupling_xolotlTestApp;

template <>
InputParameters validParams<coupling_xolotlTestApp>();

class coupling_xolotlTestApp : public MooseApp
{
public:
  coupling_xolotlTestApp(InputParameters parameters);
  virtual ~coupling_xolotlTestApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
  static void registerExecFlags(Factory & factory);
};

#endif /* COUPLING_XOLOTLTESTAPP_H */
