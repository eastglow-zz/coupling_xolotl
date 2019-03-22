coupling_xolotl
=====

"Fork coupling_xolotl" to create a new MOOSE-based application.

To use this application, MOOSE framework (https://mooseframework.inl.gov/getting_started/index.html) and Xolotl (look for URLs in following sections) should be already installed. When you install PETSc for MOOSE on Linux machine, petsc-3.10.2 (https://github.com/petsc/petsc , commit#: 2960cdc8ad2c61e9f0df8570deb766ef301642a0) which is the same one for Xolotl is recommended.

Xolotl is built into a dynamic library, and each branch of this application uses different method to implement the classes and functions from Xolotl.

'master' branch uses dynamic linking method.

To use the dynamic linking method, please install Xolotl using 'mooseApps-coupling' branch of https://github.com/ORNL-Fusion/xolotl.git.

'externCinterface' branch uses dynamic loading method, and was tested on both MacOS and Linux (HiPerGator2, University of Florida).

To use the dynamic loading method, please install Xolotl using 'mooseApps-coupling' branch of https://github.com/eastglow-zz/xolotl.git.

**You can install multiple versions of Xolotl on the same machine.

For the minimum effort to use this application, installation path of Xolotl is recommended at "../xolotl-build".

(For the record) Xolotl configuration command line prompt used:

CXX=/opt/moose/mpich-3.2/clang-6.0.1/bin/mpicxx CC=/opt/moose/mpich-3.2/clang-6.0.1/bin/mpicc FC="" PETSC_DIR=/Users/donguk.kim/myinstall/petscworkdir/petsc-clang HDF5_ROOT=/opt/hdf5 BOOST_ROOT=/opt/boost cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=no -DBUILD_SHARED_LIBS=yes  -DCMAKE_POSITION_INDEPENDENT_CODE=True -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=true /Users/donguk.kim/myinstall/xolotl-dylib/xolotl/xolotl

(Note on PETSc versions for externCinterface branch)

On MacOS:

For MOOSE: petsc-3.9.4

For Xolotl: petsc-3.10.2 (https://github.com/petsc/petsc , commit#: 2960cdc8ad2c61e9f0df8570deb766ef301642a0)

On Linux(HiPerGator2):

For both MOOSE and Xolotl:

petsc-3.10.2 (https://github.com/petsc/petsc , commit#: 2960cdc8ad2c61e9f0df8570deb766ef301642a0) 
