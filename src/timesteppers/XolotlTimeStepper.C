#include <XolotlTimeStepper.h>
#include <interface.h>
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
//    _func_name(getParam<std::string>("function_name"))
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

    XolotlInterface interface;

    interface.printSomething();

    int argc = 3;
    char ** argv = new char*[argc];
    std::string parameterFile = "crap";
    argv[0] = new char[parameterFile.length() + 1];
    strcpy(argv[0], parameterFile.c_str());
    //parameterFile = "/Users/sophie/MOOSE/moose/coupling_nonmoose_test/param.txt";
    parameterFile = "/Users/donguk.kim/projects/coupling_xolotl/params_NE.txt";
    argv[1] = new char[parameterFile.length() + 1];
    strcpy(argv[1], parameterFile.c_str());
    argv[2] = 0; // null-terminate the array

    auto solver = interface.initializeXolotl(argc,
                                             argv, MPI_COMM_WORLD);
    interface.solveXolotl(solver);
    interface.finalizeXolotl(solver);


  ////////////////
  /*
  printf("Executing internal custum routine\n");
  for (unsigned int i = 1; i <= 10; i++)
  {
    printf("%d\n",i);
  }
  */
  ////////////////
//  using std::cout;
//  using std::cerr;
//  printf("Calling an external function\n");
//  //void* handle = dlopen("/Users/donguk.kim/projects/coupling_nonmoose_test/external_app/static/lib/wrapper.so", RTLD_LAZY);
//  void* handle = dlopen(_ext_lib_path_name.c_str(), RTLD_LAZY);
//  //the library path above alse can be provided by the input file
//
//  if (!handle) {
//    cerr << "Cannot open library: " << dlerror() << '\n';
//  }
//
//  typedef void (*importfunc_t)();
//
//  dlerror();
//  importfunc_t func_loaded = (importfunc_t) dlsym(handle, _func_name.c_str());
//  const char *dlsym_error = dlerror();
//  if (dlsym_error) {
//    cerr << "Cannot load symbol " << _func_name.c_str() << ": " << dlsym_error << '\n';
//    dlclose(handle);
//  }else{
//    func_loaded();
//
//    cout << "Closing library...\n";
//    dlclose(handle);
//  }

  ////////////////
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
