#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <Fortran compiler>" >&2
    exit 2
fi

compiler=$1
if ! command -v fobis >/dev/null 2>&1; then
    echo "FoBiS.py is required to build VTKFortran." >&2
    echo "Install it as described in external/VTKFortran/README.md:" >&2
    echo "pip install FoBiS.py" >&2
    exit 1
fi

vtk_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/VTKFortran" && pwd)"
cd "$vtk_dir"

fobis fetch
cmake -S . -B build \
    -DBUILD_TESTING=ON \
    -DCMAKE_Fortran_COMPILER="$compiler"
cmake --build build
ctest --test-dir build
