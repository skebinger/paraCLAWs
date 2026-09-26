# ======================== MAKE ======================== #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

# ======================================================
# paraCLAWs User Makefile
# ======================================================

# Path to the distributed solver sources and to the user's working tree.
# Set BASE_DIR to the directory containing src/ and mk/.
BASE_DIR := $(CURDIR)
USER_DIR := $(CURDIR)

################################################################
# Install location
################################################################

inst_path = $(HOME)/.local/bin

################################################################
# Compiler
################################################################

FC = mpiifx

################################################################
# MPI ranks
################################################################

num_ranks = 12

################################################################
# Optional features
################################################################

use_vtk = 0

################################################################
# Libraries
################################################################

ifeq ($(strip $(FC)),mpifort)
LAPACK_LIB = -llapack -lblas
else ifeq ($(strip $(FC)),mpiifx)
LAPACK_LIB = -qmkl
endif

MPI_DCL_LIB_gf = -L$(HOME)/lib/mpidcl_gnu/lib64 -lmpidcl
MPI_DCL_LIB_intel = -L$(HOME)/lib/mpidcl_intel/lib64 -lmpidcl

ifeq ($(use_vtk),1)
VTKFortran_gf = -L$(HOME)/lib/VTKFortran_gnu/lib -lVTKFortran -lz
VTKFortran_intel = -L$(HOME)/lib/VTKFortran_intel/lib -lVTKFortran -lz
else
VTKFortran_gf =
VTKFortran_intel =
endif

MPI_DCL_MOD_gf = -I$(HOME)/lib/mpidcl_gnu/include
MPI_DCL_MOD_intel = -I$(HOME)/lib/mpidcl_intel/include

ifeq ($(use_vtk),1)
VTKFortran_MOD_gf = -I$(HOME)/lib/VTKFortran_gnu/include
VTKFortran_MOD_intel = -I$(HOME)/lib/VTKFortran_intel/include
else
VTKFortran_MOD_gf =
VTKFortran_MOD_intel =
endif

ifeq ($(use_vtk),1)
VTKFortran_MOD_gf += -DUSE_VTK
VTKFortran_MOD_intel += -DUSE_VTK
endif

################################################################
# Compiler Flags
################################################################

####################################################################################
# GFORTRAN
####################################################################################

ifeq ($(strip $(FC)),mpifort)

FLFLAGS = -fopenmp

FCFLAGS = \
$(MPI_DCL_MOD_gf) \
$(VTKFortran_MOD_gf) \
-J$(BUILD_DIR)/mod \
-c \
-cpp \
-ffree-line-length-none \
-fopenmp \
-O3 \
-funroll-loops \
-ftree-vectorize \
-flto=3 \
-march=native \
-finline-functions \
-fimplicit-none

FDFLAGS = -Wall -Wextra -Wconversion -pedantic

FLIBFLAGS = \
$(LAPACK_LIB) \
$(MPI_DCL_LIB_gf) \
$(VTKFortran_gf)

endif

####################################################################################
# IFX
####################################################################################

ifeq ($(strip $(FC)),mpiifx)

FLFLAGS = -qopenmp

FCFLAGS = \
$(MPI_DCL_MOD_intel) \
$(VTKFortran_MOD_intel) \
-module $(BUILD_DIR)/mod \
-c \
-cpp \
-O3 \
-funroll-loops \
-ftz \
-qopenmp \
-xHost \
-implicitnone \
-fp-model=precise

FDFLAGS = -traceback

FLIBFLAGS = \
$(LAPACK_LIB) \
$(MPI_DCL_LIB_intel) \
$(VTKFortran_intel)

endif

################################################################
# Program Name
################################################################

PROGRAM = paraCLAWs

################################################################
# Include Shared Build Rules
################################################################

include $(BASE_DIR)/mk/rules.mk