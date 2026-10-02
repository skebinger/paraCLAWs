#!/usr/bin/env python3

# ======================= python ======================= #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

"""Control settings for solver input generation.

This module defines the fields written to `control.in`.
The values are written in the order given by `FIELD_ORDER`.
"""

from dataclasses import dataclass


@dataclass
class ControlSettings:

    ###################################################################################
    # General settings
    ###################################################################################
    # number of equations solved
    number_of_equations: int = 3

    # number of dimensions; CURRENTLY NOT USED
    number_of_dimensions: int = 2 

    # number of ghost cells
    number_of_ghost_cells: int = 2

    # number of auxiliary variables; determines the size of the aux array
    number_of_auxiliary_vars: int = 10

    # number of waves in the Riemann problem; (only for hyperbolic solvers)
    number_of_waves: int = 2

    ###################################################################################
    # Run time and output settings
    ###################################################################################
    # starting time;
    time_start: float = 0 # [time_start]=s

    # final simulation time;
    time_final: float = 0.5 # [t_final]=s

    # number of snapshots to be written to disk (equal spaced)
    number_of_snapshots: int = 5

    ###################################################################################
    # Time step control settings
    ###################################################################################
    # size of the inital time step; this is also the (constant) time step used if adaptive time stepping is turned off
    initial_timestep: float = 1e-10 # [initial_timestep]=s

    # use adaptive time stepping based on a CFL condition
    adaptive_timestepping: bool = True  # True or False; False=off, True=on

    # the CFL limit for the adaptive time stepping
    CFL_limit: float = 0.1

    # minimum and maximum allowed time step sizes for the adaptive time stepping
    min_allowed_dt: float = 1e-10 # [min_allowed_dt]=s
    max_allowed_dt: float = 5e-5 # [max_allowed_dt]=s

    # maximum allowed growth and shrinkage ratios for the adaptive time stepping
    dt_max_growth_ratio: float = 1.2 # maximum allowed ratio with which the time step width is allowed to grow
    dt_min_shrink_ratio: float = 0.5 # sets the lower limit of how much the time step width is sllowd to decrease per time step update
    # Note: is the shrink ratio is exceeded, the solver will reject the time step and retry with a smaller dt

    ###################################################################################
    # Higher order flux correction
    ###################################################################################
    # Use a higher order flux correction. Currently only second order is implemented
    add_higher_order_flux_correction: bool = False  # True or False; False=off, True=on

    # Setting for limiter
    use_limtiter: bool = True  # True or False; False=off, True=on

    # Choose the limiter function; 1=minmod, 2=Van Leer, 3=superbee, 4=MC
    limiter_function: int = 1


    ###################################################################################
    # Domain settings
    ###################################################################################
    # Mesh resolution and domain dimensions
    m_xi: int = 21  # number of cells in xi direction
    m_eta: int = 21  # number of cells in eta direction

    # Size of the domain in xi and eta in the unit of choosing
    xi_min: float = -1 # lower boundary of xi
    xi_max: float = 1  # upper boundary of xi
    eta_min: float = -1  # lower boundary of eta
    eta_max: float = 1  # upper boundary of eta

    # boundary conditions
    bc_xi_lower: int = 2
    bc_xi_upper: int = 2
    bc_eta_lower: int = 2
    bc_eta_upper: int = 2

    ###################################################################################
    FIELD_ORDER = [
        "number_of_equations",
        "number_of_dimensions",
        "number_of_ghost_cells",
        "number_of_auxiliary_vars",
        "number_of_waves",
        "time_start",
        "time_final",
        "number_of_snapshots",
        "initial_timestep",
        "adaptive_timestepping",
        "CFL_limit",
        "min_allowed_dt",
        "max_allowed_dt",
        "dt_max_growth_ratio",
        "dt_min_shrink_ratio",
        "add_higher_order_flux_correction",
        "use_limtiter",
        "limiter_function",
        "m_xi",
        "m_eta",
        "xi_min",
        "xi_max",
        "eta_min",
        "eta_max",
        "bc_xi_lower",
        "bc_xi_upper",
        "bc_eta_lower",
        "bc_eta_upper"
    ]

    def values(self):
        return [getattr(self, field) for field in self.FIELD_ORDER]