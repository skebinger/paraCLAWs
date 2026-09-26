! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_solver_parameters
    implicit none

    !base solver settings
    integer, protected :: num_equations
    integer, protected :: num_waves
    integer, protected :: num_ghost
    integer, protected :: num_aux

    logical, protected :: integrate_source=.true.
    logical, protected :: add_higher_order_flux_correction=.false.
    logical, protected :: execute_post=.true.
    logical, protected :: use_limiters = .true.
    integer, protected :: limiter_method = 1

    double precision :: grav

contains

    subroutine set_solver_parameters(num_equations_,num_waves_,num_ghost_,num_aux_,add_higher_order_flux_correction_,use_limiters_,limiter_method_)
        integer :: num_equations_
        integer :: num_waves_
        integer :: num_ghost_
        integer :: num_aux_
        logical :: add_higher_order_flux_correction_
        logical :: use_limiters_
        integer :: limiter_method_

        num_equations   = num_equations_
        num_waves       = num_waves_
        num_ghost       = num_ghost_
        num_aux         = num_aux_
        add_higher_order_flux_correction = add_higher_order_flux_correction_
        use_limiters = use_limiters_
        limiter_method = limiter_method_
    end subroutine

end module MOD_solver_parameters
