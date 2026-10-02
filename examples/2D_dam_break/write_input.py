#!/usr/bin/env python3

# ======================= python ======================= #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

"""Generate input data for hyperCLAWs.

This writes the `input.data` file in the exact order expected by `setup_user_input`.
"""

from pathlib import Path
import argparse
import textwrap

from custom_settings import CustomSettings
from control_settings import ControlSettings


def format_val(val):
    if isinstance(val, bool):
        return ".true." if val else ".false."
    return val


def write_ascii_file(path: Path, values):
    with path.open("w", encoding="utf-8") as f:
        for val in values:
            f.write(f"{format_val(val)}\n")


def write_input_files(custom, control, custom_path: Path, control_path: Path):
    write_ascii_file(custom_path, custom.values())
    write_ascii_file(control_path, control.values())


def describe():
    return textwrap.dedent(
        """
        custom.in fields:

        control.in fields:
          number_of_equations        number of equations solved
          number_of_dimensions       number of dimensions; CURRENTLY NOT USED
          number_of_ghost_cells      number of ghost cells
          number_of_auxiliary_vars   number of auxiliary variables; determines the size of the aux array
          number_of_waves            number of waves in the Riemann problem; (only for hyperbolic solvers)
          time_start                 initial simulation time [time_start]=s
          time_final                 final simulation time [time_final]=s
          number_of_snapshots        number of snapshots to be written to disk
          initial_timestep           initial/constant time step [initial_timestep]=s
          adaptive_timestepping      use adaptive time stepping (True/False)
          CFL_limit                  CFL limit for adaptive time stepping
          min_allowed_dt             minimum allowed time step [min_allowed_dt]=s
          max_allowed_dt             maximum allowed time step [max_allowed_dt]=s
          dt_max_growth_ratio        maximum allowed growth ratio for adaptive time stepping
          dt_min_shrink_ratio        lower limit for allowed shrinkage ratio for adaptive time stepping
          add_higher_order_flux_correction  use higher order flux correction (True/False)
          use_limiter                use a limiter for the higher order flux correction (True/False)
          limiter_function           1=minmod, 2=Van Leer, 3=superbee, 4=MC
          m_xi                       number of cells in xi direction
          m_eta                      number of cells in eta direction
          xi_min                     lower xi boundary
          xi_max                     upper xi boundary
          eta_min                    lower eta boundary
          eta_max                    upper eta boundary
          bc_xi_lower                boundary condition in xi lower boundary
          bc_xi_upper                boundary condition in xi upper boundary
          bc_eta_lower               boundary condition in eta lower boundary
          bc_eta_upper               boundary condition in eta upper boundary
        """
    )


def print_defaults(custom, control):
    for name, value in custom.__dict__.items():
        print(f"{name} = {value}")
    for name, value in control.__dict__.items():
        print(f"{name} = {value}")


def main():
    parser = argparse.ArgumentParser(
        description="Create Fortran solver custom.in and control.in from separate parameter modules."
    )
    parser.add_argument(
        "--custom-output",
        default="custom.in",
        help="Output path for custom settings",
    )
    parser.add_argument(
        "--control-output",
        default="control.in",
        help="Output path for solver control settings",
    )
    parser.add_argument(
        "--describe",
        action="store_true",
        help="Show user-friendly input descriptions",
    )
    parser.add_argument("--defaults", action="store_true",
                        help="Show default values")
    args = parser.parse_args()

    custom = CustomSettings()
    control = ControlSettings()

    if args.describe:
        print(describe())
        return

    if args.defaults:
        print_defaults(custom, control)
        return

    write_input_files(custom, control, Path(
        args.custom_output), Path(args.control_output))
    print(
        f"Wrote {args.custom_output} and {args.control_output} with "
        f"{len(custom.values())} and {len(control.values())} values."
    )


if __name__ == "__main__":
    main()
