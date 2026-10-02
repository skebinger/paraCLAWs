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

The default `Makefile` configures the compiler via `FC`, sets `num_ranks` for MPI runs. The executable and input writer can also be installed to a user-defined directory such as `~/.local/bin` with `make install`.

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

```bash
make all
```

Useful make targets are:

```bash
make help
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

The build system supports a base source tree plus optional user-defined Fortran sources. `USER_DIR` can point at a separate working directory, and the generated source list combines the base solver files with any `.f90` files found there. If a user file declares the same Fortran module or submodule as a base source, it replaces that base source in the final build; otherwise it is added as an additional source.

This is useful for custom equations, source terms, initial conditions, or boundary logic without copying the entire base solver tree into a user directory.

## Notes

- Output is written into `output/` during execution.
- The project uses `control.in` and `custom.in` as the runtime input files for the solver configuration.
- `write_input.py --describe` prints a user-friendly summary of the expected input fields.
