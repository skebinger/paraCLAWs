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
# Set BASE_DIR to the directory containing src/ and mk/. (propably the location where you cloned the paraCLAWs repository)
BASE_DIR := ~/lib/git/paraCLAWs
# Set USER_DIR to the directory containing the user's working tree. (might be the same as BASE_DIR or a different location where you want to build the code)
USER_DIR := ~/lib/git/paraCLAWs/examples/2D_dam_break/cartesian

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

num_ranks = 5

################################################################
# Optional features
################################################################

use_vtk = 1

################################################################
# Libraries
################################################################

ifeq ($(strip $(FC)),mpifort)
LAPACK_LIB = -llapack -lblas
else ifeq ($(strip $(FC)),mpiifx)
LAPACK_LIB = -qmkl
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
$(MPI_DCL_MOD) \
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
$(MPI_DCL_LIB) \
$(VTKFortran_gf)

endif

####################################################################################
# IFX
####################################################################################

ifeq ($(strip $(FC)),mpiifx)

FLFLAGS = -qopenmp

FCFLAGS = \
$(MPI_DCL_MOD) \
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
$(MPI_DCL_LIB) \
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
