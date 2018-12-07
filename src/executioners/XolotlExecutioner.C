
#include "XolotlExecutioner.h"
#include <XolotlDLinterface.h>  // Xolotl interface
#include <dlfcn.h>

#ifdef __APPLE__
#define ISSTANDALONE true
#elif __linux__
#define ISSTANDALONE false
#endif

template<>
InputParameters validParams<XolotlExecutioner>()
{
  InputParameters params = validParams<Transient>();
  params.addParam<std::string>("library_path_name",
                               "default",
                               "Name with the path for the dynamic library to load");
  params.addParam<std::string>("XolotlInput_path_name",
                               "default",
                               "Name with the path for the Xolotl input file");
  return params;
}

XolotlExecutioner::XolotlExecutioner(const InputParameters & parameters) :
    Transient(parameters),
    _ext_lib_path_name(getParam<std::string>("library_path_name")),
    _xolotl_input_path_name(getParam<std::string>("XolotlInput_path_name"))
{

}

// This method is called only one time at the start of the entire simulation.
void XolotlExecutioner::init() {
  Transient::init();

  // Library loading
  using std::cout;
  using std::cerr;
  void* handle = dlopen(_ext_lib_path_name.c_str(), RTLD_LAZY);
  if (!handle) {
    cerr << "Cannot open library: " << dlerror() << '\n';
  }
  dlerror();

  // Class constructer & destructor loading
  typedef XolotlDLinterface* create_t();
  typedef void destroy_t(XolotlDLinterface*);

  create_t* create_interface = (create_t*) dlsym(handle, "create");
  const char* dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol create: " << dlsym_error << '\n';
  }

  destroy_t* destroy_interface = (destroy_t*) dlsym(handle, "destroy");
  dlsym_error = dlerror();
  if (dlsym_error) {
      cerr << "Cannot load symbol destroy: " << dlsym_error << '\n';
  }

  // create an instance of the class
  XolotlDLinterface* interface = create_interface();

  int argc = 3;
  char ** argv = new char*[argc];
  std::string parameterFile = "crap";
  argv[0] = new char[parameterFile.length() + 1];
  strcpy(argv[0], parameterFile.c_str());
  argv[1] = new char[_xolotl_input_path_name.length() + 1];
  strcpy(argv[1], _xolotl_input_path_name.c_str());
  argv[2] = 0; // null-terminate the array

  // MPI_Init(&argc, &argv);
  // Initialize it
  //auto solver = interface->initializeXolotl(argc, argv, MPI::COMM_WORLD);
  //auto solver = interface->initializeXolotl(argc, argv, MPI_COMM_WORLD);
  auto solver = interface->initializeXolotl(argc, argv, MPI_COMM_WORLD, ISSTANDALONE); // 'isStandalone=false' doesn't work on MacOS. But, 'true' works. I don't know why.
  // std::shared_ptr<xolotlSolver::PetscSolver> solver = interface->initializeXolotl(argc, argv, MPI_COMM_WORLD);
  // Run the solve
  interface->solveXolotl(solver);
  // Finalize the run
  // interface->finalizeXolotl(solver);
  // interface->finalizeXolotl(solver, ISSTANDALONE); // 'isStandalone=false' doesn't work on MacOS. But, 'true' works. I don't know why.
  interface->finalizeXolotl(solver, false); // 'isStandalone=false' doesn't work on MacOS. But, 'true' works. I don't know why.

  // destroy the class
  destroy_interface(interface);

  // unload the triangle library
  dlclose(handle);

  // MPI_Finalize();

  ////////////////
}
