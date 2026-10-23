#!/usr/bin/env bash

# ======================== bash ======================== #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #


# ----------------------------------------------------------------------
# Opens ParaView to visualize the output of the solver.
# ----------------------------------------------------------------------

set -euo pipefail

output_dir="$PWD/output"

if [[ ! -d "$output_dir" ]]; then
    printf 'Solver output directory not found: %s\nRun the solver first with VTK output enabled.\n' "$output_dir" >&2
    exit 1
fi

pvd_files=("$output_dir"/*.pvd)
if [[ ! -e "${pvd_files[0]}" ]]; then
    printf 'No ParaView collection (*.pvd) found in %s\nRun the solver with VTK output enabled.\n' "$output_dir" >&2
    exit 1
fi

if [[ ${#pvd_files[@]} -gt 1 ]]; then
    printf 'Multiple ParaView collections found in %s; expected exactly one:\n' "$output_dir" >&2
    printf '  %s\n' "${pvd_files[@]}" >&2
    exit 1
fi

if ! command -v paraview >/dev/null 2>&1; then
    printf 'ParaView executable not found on PATH.\n' >&2
    exit 127
fi

exec paraview "${pvd_files[0]}"
