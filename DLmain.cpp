#include <XolotlDLinterface.h>
#include <iostream>
#include <dlfcn.h>
#include <petscsys.h>

double* build_xolotl_axis(int nsize, double dl)
{
  double *tmpbuff;
  tmpbuff = new double[nsize];
  for (int i = 0; i < nsize; i++) {
    tmpbuff[i] = (double) i * dl;
  }
  return tmpbuff;
}

double* allocate_buffer(int nsize)
{
	double *buff;
	buff = new double[nsize];
	return buff;
}

void init_buffer(double *buff, int nsize, double val)
{
	for (int i = 0; i < nsize; i++){
		buff[i] = val;
	}
}

void local_index_info(int *xsreturn, int *xmreturn, XolotlDLinterface* interface, std::shared_ptr<xolotlSolver::PetscSolver> solver)
{
	int xs, xm, d;
	interface->getLocalCoordinates(solver, &xs, &xm, &d, &d, &d, &d, &d, &d, &d);

	*xsreturn = xs;
	*xmreturn = xm;
}

void reinit_with_local_data(double *buff, int GlobSize, int xs, int xm, XolotlDLinterface* interface, std::shared_ptr<xolotlSolver::PetscSolver> solver)
{
	std::vector<std::vector<std::vector<double> > > *tmp = interface->getLocalXeConc(solver);
	for (int i = 0; i < GlobSize; i++){
		buff[i] = 0.0;
	}
	for (int i = xs; i < xs + xm; i++) {
		buff[i] = tmp->at(i-xs)[0][0];
	}
}

// void merge_global_data(double *buff, int GlobSize, MPI_Comm comm)
// {
// 	for (int i = 0; i < GlobSize; i++) {
// 		double tmpsum = 0.0;
// 		MPI_Allreduce(&buff[i], &tmpsum, 1, MPI_DOUBLE, MPI_SUM, comm);
// 		buff[i] = tmpsum;
// 	}
// }

void fileout(char *fname, double *x, double *data, int xs, int xf)
{
	FILE *foutput;
	foutput = fopen(fname, "w");
	for (int i = xs; i < xf; i++) {
		fprintf(foutput, "%lf %lf\n", x[i], data[i]);
	}
	fclose(foutput);
}

//! Main program
int main(int argc, char **argv) {
	using std::cout;
	using std::cerr;

	// Loading libxolotlInter.so
	void* handle = dlopen("/Users/donguk.kim/projects/xolotl-build/lib/libxolotlInter.dylib",RTLD_LAZY);
	if (!handle) {
		cerr << "cannot load library: " << dlerror() << '\n';
		return 1;
	}
	dlerror();

	// Loading XolotlInterface class
	//XIallocator_t* interface = (XIallocator_t*) dlsym(interfacelib, "XIallocate");
	typedef XolotlDLinterface* create_t();
	typedef void destroy_t(XolotlDLinterface*);
        create_t* create_interface = (create_t*) dlsym(handle, "create");
	const char* dlsym_error = dlerror();
	if (dlsym_error) {
		cerr << "cannot load symbol create: " << dlsym_error << '\n';
		return 1;
	}

	destroy_t* destroy_interface = (destroy_t*) dlsym(handle, "destroy");
	dlsym_error = dlerror();
	if (dlsym_error) {
		cerr << "cannot load symbol destroy: " << dlsym_error << '\n';
		return 1;
	}

  // Creating an instance of the class
	XolotlDLinterface* interface = create_interface();

	// Initialize MPI
	MPI_Init(&argc, &argv);
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// Initialize it
	auto solver = interface->initializeXolotl(argc, argv, MPI_COMM_WORLD, true);

	int nx, dummy;
	double dx, lx, ldummy;
	interface->getXolotlGlobalGridInfo(&dummy, &nx, &dummy, &dummy, &dx, &ldummy, &ldummy, &lx, &ldummy, &ldummy, argc, argv);
  double *x, *conc;
	x = build_xolotl_axis(nx, dx);
	conc = allocate_buffer(nx);
	init_buffer(conc, nx, 0);
  int xs, xm;
	local_index_info(&xs, &xm, interface, solver);
	printf("myrank = %d, nx = %d, xs = %d, xm = %d\n",myrank, nx, xs, xm);

	// Run the solve
	interface->solveXolotl(solver);
	reinit_with_local_data(conc, nx, xs, xm, interface, solver);
	// merge_global_data(conc, nx, MPI_COMM_WORLD);

	// for (int i = 0; i < nx; i++) {
	// 	double tmpsum = 0.0;
	// 	MPI_Allreduce(&conc[i], &tmpsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// 	conc[i] = tmpsum;
	// }

	char filename[100];
	// if (myrank == 1) {
	// 	printf("sizeof x[] = %lu\n", sizeof(x)/sizeof(double));
	// 	printf("sizeof conc[] = %lu\n", sizeof(conc)/sizeof(double));
	// 	sprintf(filename, "./dataout_%02d.dat", myrank);
	// 	fileout(filename, x, conc, xs, xs+xm);
	// }
	sprintf(filename, "./dataout_%02d.dat", myrank);
	fileout(filename, x, conc, xs, xs+xm);

	// Finalize the run
	interface->finalizeXolotl(solver);

	delete[] x;
	delete[] conc;

	// Finalize MPI
	MPI_Finalize();


	// // Initialize MPI
	// MPI_Init(&argc, &argv);
	//
	// // Create an interface to control the solver
	// XolotlInterface interface;
	//
	// // Initialize it
	// auto solver = interface.initializeXolotl(argc, argv);
	// // Run the solve
	// interface.solveXolotl(solver);
	// // Finalize the run
	// interface.finalizeXolotl(solver);
	//
	// // Finalize MPI
	// MPI_Finalize();

	return EXIT_SUCCESS;
}
