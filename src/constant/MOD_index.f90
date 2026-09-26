! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_index
    !! Module containing index definitions.
    !! For now only contains the aux. array.
    implicit none

    ! auxiliary array
    integer, parameter :: idx_aux_uw = 1
    integer, parameter :: idx_aux_vw = 2
    integer, parameter :: idx_aux_ui = 3
    integer, parameter :: idx_aux_vi = 4
    integer, parameter :: idx_aux_ut = 5
    integer, parameter :: idx_aux_vt = 6
    integer, parameter :: idx_aux_du = 7
    integer, parameter :: idx_aux_dv = 8
    integer, parameter :: idx_aux_tau_x = 9
    integer, parameter :: idx_aux_tau_y = 10

contains

end module MOD_index
