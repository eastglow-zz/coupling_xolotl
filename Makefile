###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES         := no

CHEMICAL_REACTIONS  := no
CONTACT             := no
FLUID_PROPERTIES    := no
HEAT_CONDUCTION     := no
MISC                := no
NAVIER_STOKES       := no
PHASE_FIELD         := yes
RDG                 := no
RICHARDS            := no
SOLID_MECHANICS     := no
STOCHASTIC_TOOLS    := no
TENSOR_MECHANICS    := no
XFEM                := no
POROUS_FLOW         := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################
#External non-MOOSE shared library linking - not necessary
extlibpath          :=$(XOLOTL_PATH)/lib
extheaderpath       :=$(XOLOTL_PATH)/include
ADDITIONAL_INCLUDES:= -I$(extheaderpath)
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := coupling_xolotl
BUILD_EXEC         := yes
GEN_REVISION       := no
include            $(FRAMEWORK_DIR)/app.mk
export EXTERNAL_FLAGS += -std=c++11 -L$(extlibpath) -lxolotlInter

###############################################################################
# Additional special case targets should be added here
xo:
	@mkdir xolotl/build
	@cd xolotl/build; PETSC_DIR=/home/sophie/MOOSE/projects/moose/petsc/ PETSC_ARCH=arch-moose \
	CXX=mpicxx CC=mpicc HDF5_ROOT=/home/sophie/Code/hdf5-install/ \
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx \
	-DBUILD_SHARED_LIBS=yes -DCMAKE_CXX_FLAGS_RELEASE="-o3 -fPIC" ../xolotl; \
	make; make install
