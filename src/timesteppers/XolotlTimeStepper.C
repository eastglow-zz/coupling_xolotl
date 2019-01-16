#include <XolotlTimeStepper.h>
#include <dlfcn.h>

template<>
InputParameters validParams<XolotlTimeStepper>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addParam<Real>("dt", 1.0, "Size of the time step");
  params.addParam<std::string>("library_path_name",
                               "default",
                               "Name with the path for the dynamic library to load");
//  params.addParam<std::string>("function_name",
//                               "default",
//                               "Name of the function to load");
  return params;
}

XolotlTimeStepper::XolotlTimeStepper(const InputParameters & parameters) :
    TimeStepper(parameters),
    _dt(getParam<Real>("dt")),
    _ext_lib_path_name(getParam<std::string>("library_path_name"))
{
    // Initialize the step number
    stepNumber = 0;

    int argc = 3;
    char ** argv = new char*[argc];
    std::string parameterFile = "crap";
    argv[0] = new char[parameterFile.length() + 1];
    strcpy(argv[0], parameterFile.c_str());
//    parameterFile = "/home2/bqo/MOOSE/projects/coupling_xolotl_std/params_NE.txt";
    parameterFile = "/Users/donguk.kim/projects/coupling_xolotl/params_NE_3D.txt";    
    argv[1] = new char[parameterFile.length() + 1];
    strcpy(argv[1], parameterFile.c_str());
    argv[2] = 0; // null-terminate the array

    solver = interface.initializeXolotl(argc,
                                                                 argv, MPI_COMM_WORLD, true);
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
    // Increment the step number
    stepNumber++;

    // Set the time we want to reach
    interface.setTimes(solver, 5.0e4 * stepNumber, 5.0e3);

    // Run the solver
    interface.solveXolotl(solver);
    
    // Get the local rate
    auto localRate = interface.getLocalXeRate(solver);
    // Print its dimensions
    std::cout << "The vector is: " << localRate->size() << " in the X direction, " << localRate->at(0).size()
	<< " in the Y direction, and " << localRate->at(0)[0].size() << " in the Z direction." << std::endl; 
}

void
XolotlTimeStepper::postExecute()
{
    interface.finalizeXolotl(solver, false);

    solver.reset();
}

/* Indication of whether the Monte Carlo solve converged - assume this will
   always be true. */
bool
XolotlTimeStepper::converged()
{
  return true;
}
