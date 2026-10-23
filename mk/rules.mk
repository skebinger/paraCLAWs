# ======================== MAKE ======================== #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

# ======================== Make ======================== #
# Shared Build Rules
# ======================================================

.DEFAULT_GOAL := help

override BASE_DIR := $(abspath $(patsubst ~/%,$(HOME)/%,$(BASE_DIR)))
override USER_DIR := $(abspath $(patsubst ~/%,$(HOME)/%,$(USER_DIR)))

MPIDCL_INSTALL_DIR := $(BASE_DIR)/external/mpidcl/install
MPI_DCL_MOD := -I$(MPIDCL_INSTALL_DIR)/include
MPI_DCL_LIB := -L$(MPIDCL_INSTALL_DIR)/lib64 -lmpidcl

VTKFORTRAN_BUILD_DIR := $(BASE_DIR)/external/VTKFortran/build
VTKFORTRAN_LIBS = \
$(VTKFORTRAN_BUILD_DIR)/lib64/libVTKFortran.a \
$(VTKFORTRAN_BUILD_DIR)/src/third_party/FoXy/lib64/libFoXy.a \
$(VTKFORTRAN_BUILD_DIR)/src/third_party/StringiFor/lib64/libStringiFor.a \
$(VTKFORTRAN_BUILD_DIR)/src/third_party/BeFoR64/lib64/libBeFoR64.a \
$(VTKFORTRAN_BUILD_DIR)/src/third_party/PENF/lib64/libPENF.a \
$(VTKFORTRAN_BUILD_DIR)/src/third_party/FACE/lib64/libFACE.a
VTKFORTRAN_MODS = \
-I$(VTKFORTRAN_BUILD_DIR)/modules \
-I$(VTKFORTRAN_BUILD_DIR)/src/third_party/FoXy/modules \
-I$(VTKFORTRAN_BUILD_DIR)/src/third_party/StringiFor/modules \
-I$(VTKFORTRAN_BUILD_DIR)/src/third_party/BeFoR64/modules \
-I$(VTKFORTRAN_BUILD_DIR)/src/third_party/PENF/modules \
-I$(VTKFORTRAN_BUILD_DIR)/src/third_party/FACE/modules

ifeq ($(use_vtk),1)
VTKFortran_gf = $(VTKFORTRAN_LIBS)
VTKFortran_intel = $(VTKFORTRAN_LIBS)
VTKFortran_MOD_gf = $(VTKFORTRAN_MODS) -DUSE_VTK
VTKFortran_MOD_intel = $(VTKFORTRAN_MODS) -DUSE_VTK
else
VTKFortran_gf =
VTKFortran_intel =
VTKFortran_MOD_gf =
VTKFortran_MOD_intel =
endif

################################################################
# Build directories
################################################################

BUILD_DIR := $(USER_DIR)/build

OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
BIN_DIR := $(BUILD_DIR)/bin

TARGET := $(BIN_DIR)/$(PROGRAM)

################################################################
# Source files (ORDER MATTERS FOR FORTRAN MODULES)
################################################################
# Order is important due to module dependencies.

BASE_SRCS = \
$(BASE_DIR)/src/constant/MOD_index.f90 \
$(BASE_DIR)/src/constant/MOD_constants.f90 \
$(BASE_DIR)/src/constant/MOD_solver_parameters.f90 \
$(BASE_DIR)/src/user_subroutines/MOD_user_data.f90 \
$(BASE_DIR)/src/utility/MOD_utils.f90 \
$(BASE_DIR)/src/utility/MOD_error_handling.f90 \
$(BASE_DIR)/src/utility/MOD_look_up_table.f90 \
$(BASE_DIR)/src/utility/MOD_mathematical_functions.f90 \
$(BASE_DIR)/src/utility/SMOD_derivatives.f90 \
$(BASE_DIR)/src/utility/SMOD_eigenvalue_problem.f90 \
$(BASE_DIR)/src/utility/SMOD_matrix_operations.f90 \
$(BASE_DIR)/src/solver/MOD_work_arrays.f90 \
$(BASE_DIR)/src/time_control/MOD_time_control.f90 \
$(BASE_DIR)/src/domain/MOD_domain.f90 \
$(BASE_DIR)/src/domain/SMOD_domain_allocation.f90 \
$(BASE_DIR)/src/domain/SMOD_mesh_manipulation.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_boundary_conditions.f90 \
$(BASE_DIR)/src/user_subroutines/MOD_user_subroutines.f90 \
$(BASE_DIR)/src/solver/SMOD_initialize_work_arrays.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_user_post.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_setup.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_set_initial_conditions.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_source_function.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_coordinate_mapping.f90 \
$(BASE_DIR)/src/solver/MOD_flux_calculation.f90 \
$(BASE_DIR)/src/solver/SMOD_approximate_Riemann_solvers.f90 \
$(BASE_DIR)/src/solver/SMOD_jakobian.f90 \
$(BASE_DIR)/src/solver/MOD_flux_correction.f90 \
$(BASE_DIR)/src/solver/MOD_solve.f90 \
$(BASE_DIR)/src/utility/MOD_IO.f90 \
$(BASE_DIR)/src/main.f90

################################################################
# Combined source list
################################################################

SOURCE_LIST_MK := $(BUILD_DIR)/source_list.mk

.PHONY: FORCE
FORCE:

$(SOURCE_LIST_MK): FORCE $(BASE_DIR)/mk/generate_source_list.py
	@mkdir -p $(dir $@)
	@python3 $(BASE_DIR)/mk/generate_source_list.py \
		--user-dir "$(USER_DIR)" \
		--base-dir "$(BASE_DIR)" \
		--output "$@" \
		$(foreach source,$(BASE_SRCS),--base-source "$(source)")

include $(SOURCE_LIST_MK)

BASE_OBJS := \
$(patsubst $(BASE_DIR)/%.f90,$(OBJ_DIR)/base/%.o,$(BASE_FINAL_SRCS))

USER_OBJS := \
$(patsubst $(USER_DIR)/%.f90,$(OBJ_DIR)/user/%.o,$(USER_FINAL_SRCS))

OBJS := $(foreach source,$(FINAL_SRCS),\
$(if $(filter $(source),$(USER_FINAL_SRCS)),\
$(patsubst $(USER_DIR)/%.f90,$(OBJ_DIR)/user/%.o,$(source)),\
$(patsubst $(BASE_DIR)/%.f90,$(OBJ_DIR)/base/%.o,$(source))))

################################################################
# Build targets
################################################################

all: directories $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FLFLAGS) $(FDFLAGS) -o $@ $^ $(FLIBFLAGS)

new:
	$(MAKE) clean
	$(MAKE) dep
	$(MAKE) all

################################################################
# Compilation rules
################################################################

$(OBJ_DIR)/base/%.o: $(BASE_DIR)/%.f90 | directories
	@mkdir -p $(dir $@)
	$(FC) $(FCFLAGS) $(FDFLAGS) -c $< -o $@

$(OBJ_DIR)/user/%.o: $(USER_DIR)/%.f90 | directories
	@mkdir -p $(dir $@)
	$(FC) $(FCFLAGS) $(FDFLAGS) -c $< -o $@

################################################################
# Directory setup
################################################################

DIRS := $(BUILD_DIR) $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)

directories:
	@mkdir -p $(DIRS)

################################################################
# Utilities
################################################################

.PHONY: dep
dep:
	@case "$(strip $(FC))" in \
		mpifort) mpi_flavor=openmpi ;; \
		mpiifx) mpi_flavor=intelmpi ;; \
		*) echo "Unsupported Fortran compiler '$(FC)' (expected mpifort or mpiifx)"; exit 1 ;; \
	esac; \
	cd "$(BASE_DIR)/external/mpidcl" && \
		./makeall.sh --type Release --shared OFF --mpi "$$mpi_flavor" \
			--build-lib ON --build-tests OFF
ifeq ($(use_vtk),1)
	@bash "$(BASE_DIR)/external/build_vtkfortran.sh" "$(FC)"
endif

run: directories $(TARGET) input
	@./set_env.sh && $(TARGET)

mpirun: directories $(TARGET) input
	@source ./set_env.sh && mpirun --bind-to socket -np $(num_ranks) $(TARGET)

create_inst_path:
	@mkdir -p $(inst_path)

install: create_inst_path $(TARGET)
	@echo "======================================================"
	@echo "Installing $(TARGET) to $(inst_path)"
	@cp $(TARGET) $(inst_path)
	@echo "Installing the input writer to $(inst_path)"
	@cp write_input.py $(inst_path)
	@echo "======================================================"

doc:
	@ford ford_doc_config.md

input:
	@python3 write_input.py

input_help:
	@python3 write_input.py --describe

clean:
	@rm -rf $(BUILD_DIR)

################################################################
# Debugging
################################################################

# If something doesn't work right, use 'make debug' to
# show what each variable contains.
debug:
	@echo "======================================================"
	@echo "BASE_DIR = $(BASE_DIR)"
	@echo "======================================================"
	@echo "USER_DIR = $(USER_DIR)"
	@echo "======================================================"
	@echo "BASE_SRCS = $(BASE_SRCS)"
	@echo "======================================================"
	@echo "BASE_FINAL_SRCS = $(BASE_FINAL_SRCS)"
	@echo "======================================================"
	@echo "USER_FINAL_SRCS = $(USER_FINAL_SRCS)"
	@echo "======================================================"
	@echo "FINAL_SRCS = $(FINAL_SRCS)"
	@echo "======================================================"
	@echo "BASE_OBJS = $(BASE_OBJS)"
	@echo "======================================================"
	@echo "USER_OBJS = $(USER_OBJS)"
	@echo "======================================================"
	@echo "OBJS = $(OBJS)"
	@echo "======================================================"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "======================================================"
	@echo "TARGET = $(TARGET)"
	@echo "======================================================"
	@echo "COMPILER FLAGS = $(FCFLAGS)"
	@echo "======================================================"
	@echo "LINKER FLAGS = $(FLFLAGS)"
	@echo "======================================================"
	@echo "DEBUG FLAGS = $(FDFLAGS)"
	@echo "======================================================"
	@echo "LIBRARY FLAGS = $(FLIBFLAGS)"
	@echo "======================================================"

print_sources:
	@printf "%s\n" $(FINAL_SRCS)

################################################################
# Help
################################################################

help:
	@echo "Usage: make [target]"
	@echo "Targets:"
	@echo "  all: compiles and links the solver"
	@echo "  new: cleans and rebuilds the solver"
	@echo "  dep: builds mpidcl; also builds VTKFortran when use_vtk=1"
	@echo "  run: compiles and runs the code (with checks for recompilation)"
	@echo "  mpirun: compiles and runs the code with MPI (with checks for recompilation)"
	@echo "  install: compiles and installs the executable to $(inst_path)"
	@echo "  doc: generates documentation using Ford"
	@echo "  input: generates a default input file using write_input.py"
	@echo "  input_help: describes the input file format using write_input.py"
	@echo "  debug: prints the active source, object, and compiler variables"
	@echo "  clean: removes the build directory and its generated files"
	@echo "  print_sources: prints the final concatenated list of source files"

.PHONY: \
all \
new \
dep \
clean \
run \
mpirun \
install \
debug \
help \
print_sources\