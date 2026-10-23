# paraCLAWs

paraCLAWs is a parallel solver for hyperbolic systems of conservation laws. The code is built around a 2D finite-volume formulation, uses MPI for distributed-memory parallelism and OpenMP for threaded execution, and includes adaptive time stepping, boundary conditions, source terms, and higher-order correction machinery.

The core solver lives in `src/`, the shared build logic lives in `mk/`, and the generated input files are produced by the Python helpers in the project root.

## Project layout

- `src/` – Fortran solver source files
- `mk/` – shared Makefile rules for building the solver
- `write_input.py` – generates solver input files from the Python settings modules
- `control_settings.py` – default control parameters writing `control.in`
- `custom_settings.py` – user/custom settings writing `custom.in`
- `set_env.sh` – environment configuration for MPI/OpenMP runtime settings (useful for HPC clusters)
- `Makefile` – project-level build entry points

## Requirements

The project expects a Fortran MPI toolchain, for example:

- `mpiifx` (Intel oneAPI) or `mpifort` (GNU)
- OpenMP support
- MPI runtime and development libraries
- BLAS/LAPACK or MKL, depending on compiler configuration
- `mpidcl` from https://github.com/skebinger/mpidcl
- optional VTK libraries from https://github.com/szaghi/VTKFortran if the corresponding features are enabled
- optional `figlet` to print a header to console (purely visual)

The default `Makefile` configures the compiler via `FC`, sets `num_ranks` for MPI runs. The executable and VTK input writer can also be installed to a user-defined directory such as `~/.local/bin` with `make install`.

## Quick start

1. Edit the compiler selection in `Makefile` if needed:

```make
FC = mpiifx
# or FC = mpifort
```

2. Generate the default input files:

```bash
make input
```

This writes `control.in` and `custom.in` in the project root.

3. Build the solver:

If compiled for the first time or recompile solver and dependencies
```bash
make new
```
Only compile solver
```bash
make all
```

Useful make targets are:

```bash
make help
make dep
make run
make mpirun
make input
make input_help
make clean
```

## Running the solver

The build output is placed under `build/bin/` and the executable is named `paraCLAWs` by default.

For a local serial run:

```bash
source ./set_env.sh
./build/bin/paraCLAWs
```

For a local parallel MPI run:

```bash
source ./set_env.sh
mpirun --bind-to socket -np 12 ./build/bin/paraCLAWs
```

## User source overrides

Keep the solver repository as the base source tree, and do your case-specific work in a separate directory. You only need to copy the Fortran files you intend to change (often files from `src/user_subroutines/`), along with any additional user source files. There is no need to copy the entire solver source tree.

For example, create a workspace and copy the project Makefile and the user hooks you want to customize:

```bash
mkdir -p ~/work/my_case
cp Makefile ~/work/my_case/
cp src/user_subroutines/SMOD_setup.f90 \
   src/user_subroutines/SMOD_set_initial_conditions.f90 \
   ~/work/my_case/
cp write_input.py control_settings.py custom_settings.py \
   ~/work/my_case/
cp set_env.sh ~/work/my_case/
```

In the copied Makefile, set `BASE_DIR` to the paraCLAWs repository, then run make from the workspace. `USER_DIR` defaults to `$(CURDIR)`, so it does not need to be set separately when make is run there.

The build scans `USER_DIR` for `.f90` files and combines them with the base solver sources. A user file that declares the same Fortran module or submodule as a base source replaces that source; otherwise, it is added to the build. Thus, copy only the hooks you need to customize, and leave the other solver sources in the base repository. The `examples/` directories demonstrate this layout.

## Notes

- Output is written into `output/` during execution.
- With `use_vtk = 1`, open `output/solution.pvd` in ParaView to browse the saved snapshots using its time controls. Each timestep's `.pvts` file references the individual rank `.vts` pieces.
- The project uses `control.in` and `custom.in` as the runtime input files for the solver configuration.
- `write_input.py --describe` prints a user-friendly summary of the expected input fields.
- `make dep` builds MPIDCL and also builds/tests VTKFortran when `use_vtk = 1` in the active Makefile.
- When `use_vtk = 1`, `make dep` requires FoBiS.py; install it with `pip install FoBiS.py` if needed.
