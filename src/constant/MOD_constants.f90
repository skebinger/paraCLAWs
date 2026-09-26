! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_constants
    implicit none

    integer, parameter :: axis_xi = 1
    integer, parameter :: axis_eta = 2

    integer, parameter :: bc_zero_gradient = 1

    double precision, parameter :: ud_precision = 1d-10

    double precision, parameter :: PI = 4.D0*DATAN(1.D0)
    
contains
    
end module MOD_constants