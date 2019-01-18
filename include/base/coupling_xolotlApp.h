//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef COUPLING_XOLOTLAPP_H
#define COUPLING_XOLOTLAPP_H

#include "MooseApp.h"

class coupling_xolotlApp;

template <>
InputParameters validParams<coupling_xolotlApp>();

class coupling_xolotlApp : public MooseApp
{
public:
  coupling_xolotlApp(InputParameters parameters);
  virtual ~coupling_xolotlApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* COUPLING_XOLOTLAPP_H */
