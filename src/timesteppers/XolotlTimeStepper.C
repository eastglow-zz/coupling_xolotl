#include "XolotlTimeStepper.h"


template<>
InputParameters validParams<XolotlTimeStepper>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addParam<Real>("dt", 1.0, "Size of the time step");
  return params;
}

XolotlTimeStepper::XolotlTimeStepper(const InputParameters & parameters) :
    TimeStepper(parameters),
    _dt(getParam<Real>("dt"))
{
}

/* This method is a pure virtual function, so it must be redefined by any
daughter class, even if the behavior doesn't change. */
Real
XolotlTimeStepper::computeInitialDT()
{
  return _dt;
}

/* This method is a pure virtual function, so it must be redefined by any
daughter class, even if the behavior doesn't change. */
Real
XolotlTimeStepper::computeDT()
{
  return _dt;
}

void
XolotlTimeStepper::step()
{
    //
    // XolotlInterface interface;
    //
    // interface.printSomething();
    //
    // int argc = 3;
    // char ** argv = new char*[argc];
    // std::string parameterFile = "crap";
    // argv[0] = new char[parameterFile.length() + 1];
    // strcpy(argv[0], parameterFile.c_str());
    // //parameterFile = "/Users/sophie/MOOSE/moose/coupling_nonmoose_test/param.txt";
    // parameterFile = "/Users/donguk.kim/projects/coupling_xolotl/params_NE.txt";
    // argv[1] = new char[parameterFile.length() + 1];
    // strcpy(argv[1], parameterFile.c_str());
    // argv[2] = 0; // null-terminate the array
    //
    // auto solver = interface.initializeXolotl(argc,
    //                                          argv, MPI_COMM_WORLD);
    // interface.solveXolotl(solver);
    // interface.finalizeXolotl(solver);
}

void
XolotlTimeStepper::postExecute()
{
}

/* Indication of whether the Monte Carlo solve converged - assume this will
   always be true. */
bool
XolotlTimeStepper::converged()
{
  return true;
}
