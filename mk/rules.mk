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
$(BASE_DIR)/src/utility/MOD_utils.f90 \
$(BASE_DIR)/src/utility/MOD_error_handling.f90 \
$(BASE_DIR)/src/utility/MOD_look_up_table.f90 \
$(BASE_DIR)/src/utility/MOD_mathematical_functions.f90 \
$(BASE_DIR)/src/solver/MOD_work_arrays.f90 \
$(BASE_DIR)/src/user_subroutines/MOD_user_data.f90 \
$(BASE_DIR)/src/time_control/MOD_time_control.f90 \
$(BASE_DIR)/src/domain/MOD_domain.f90 \
$(BASE_DIR)/src/domain/SMOD_domain_allocation.f90 \
$(BASE_DIR)/src/domain/SMOD_mesh_manipulation.f90 \
$(BASE_DIR)/src/domain/SMOD_boundary_conditions.f90 \
$(BASE_DIR)/src/user_subroutines/MOD_user_subroutines.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_user_post.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_setup.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_set_initial_conditions.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_source_function.f90 \
$(BASE_DIR)/src/user_subroutines/SMOD_coordinate_mapping.f90 \
$(BASE_DIR)/src/solver/SMOD_initialize_work_arrays.f90 \
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

new: clean all

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
clean \
run \
mpirun \
install \
debug \
help \
print_sources\