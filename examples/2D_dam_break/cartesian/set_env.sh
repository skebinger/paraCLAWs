# ======================== bash ======================== #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

# ----------------------------------------------------------------------
# Prevent BLAS libraries from oversubscribing CPU cores
# (important when running MPI + OpenMP or other threaded codes)
# ----------------------------------------------------------------------

# Force OpenBLAS to use a single thread per process
# Without this, OpenBLAS may spawn many threads per MPI rank,
# leading to severe oversubscription and poor performance.
export OPENBLAS_NUM_THREADS=1

# Same idea for GotoBLAS (older BLAS implementation)
export GOTO_NUM_THREADS=1

# Force Intel MKL to use a single thread per process
export MKL_NUM_THREADS=1

# Disable MKL’s internal dynamic thread adjustment
# Ensures thread count stays fixed and predictable
export MKL_DYNAMIC=FALSE   # optional but recommended for MPI jobs


# ----------------------------------------------------------------------
# OpenMP runtime configuration
# ----------------------------------------------------------------------

# Number of OpenMP threads per MPI rank
# Choose this so that: (MPI ranks × OMP_NUM_THREADS) ≈ total CPU cores
export OMP_NUM_THREADS=20

# Disable dynamic adjustment of OpenMP threads
# Prevents OpenMP from changing thread count at runtime,
# which can interfere with MPI load balancing and CPU pinning.
export OMP_DYNAMIC=FALSE

# Control how OpenMP threads are pinned to CPU cores (Intel runtime)
# granularity=fine  → bind threads at the core (not socket) level
# compact           → pack threads close together for cache locality
# 1,0               → one thread per core, no skipping
export KMP_AFFINITY=granularity=fine,compact,1,0

# OpenMP loop scheduling policy
# GUIDED           → decreasing chunk sizes for better load balance
# 10               → minimum chunk size of 10 iterations
# Useful when loop iteration cost is uneven
export OMP_SCHEDULE="GUIDED,10"
