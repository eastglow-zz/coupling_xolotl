coupling_xolotl
=====

This is a [MOOSE](https://mooseframework.inl.gov/getting_started/index.html) application wrapping [Xolotl](https://github.com/ORNL-Fusion/xolotl/wiki) a cluster dynamics code.

Here is how to install this application:
```
git clone https://github.com/eastglow-zz/coupling_xolotl coupling_xolotl
cd coupling_xolotl
git submodule init
git submodule update
cd moose
./scripts/update_and_rebuild_petsc.sh --download-hdf5
./scripts/update_and_rebuild_libmesh.sh
cd ..
make
```

Tests can be run through:
```
./run_tests
```

If your machine has N cores available the installation can go faster by using:
```
MOOSE_JOBS=N ./scripts/update_and_rebuild_libmesh.sh
```
for libMesh; and:
```
make -j N
```
for the coupling code.

If you want 64bit indices support simply add the `--with-64-bit-indices` option with the PETSc script.
