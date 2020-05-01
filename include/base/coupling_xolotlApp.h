//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include "MooseApp.h"
#include <interface.h>

class coupling_xolotlApp: public MooseApp {
public:
	static InputParameters validParams();

	coupling_xolotlApp(InputParameters parameters);
	virtual ~coupling_xolotlApp();

	XolotlInterface & getInterface() {
		return _interface;
	}
	TS & getXolotlTS() {
		return _interface.getTS();
	}
	static void registerApps();
	static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

	// For restart capabilities
	std::shared_ptr<Backup> backup();
	void restore(std::shared_ptr<Backup> backup, bool for_restart = false);

private:
	XolotlInterface _interface;
};
