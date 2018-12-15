coupling_xolotl
=====

"Fork coupling_xolotl" to create a new MOOSE-based application.

Xolotl is required to be built into a dynamic library.
To build Xolotl dynamic library, please use 'mooseApps-coupling' branch of https://github.com/ORNL-Fusion/xolotl.git.

'master' branch uses dynamic linking method.

'externCinterface' branch uses dynamic loading method, and was tested on MacOS.

'externCinterfaceHPG2' branch used dynamic loading method, and was tested on the HiPerGator2 cluster of University of Florida.

To use the dynamic loading method, please use 'mooseApps-coupling' branch of https://github.com/eastglow-zz/xolotl.git.

(For the record)Xolotl configure command line prompt used:

CXX=/opt/moose/mpich-3.2/clang-6.0.1/bin/mpicxx CC=/opt/moose/mpich-3.2/clang-6.0.1/bin/mpicc FC="" PETSC_DIR=/Users/donguk.kim/myinstall/petscworkdir/petsc-clang HDF5_ROOT=/opt/hdf5 BOOST_ROOT=/opt/boost cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=no -DBUILD_SHARED_LIBS=yes  -DCMAKE_POSITION_INDEPENDENT_CODE=True -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=true /Users/donguk.kim/myinstall/xolotl-dylib/xolotl/xolotl
