! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_look_up_table
    implicit none

    type lu_tab_1D
        double precision, allocatable :: table(:)
    contains
        procedure :: create_LUT
        procedure :: deallocate_LUT
        procedure :: find_LUT
    end type

    interface
        subroutine write_LUT(lut)
            import :: lu_tab_1D
            class(lu_tab_1D) :: lut
        end subroutine
    end interface

contains

    subroutine create_LUT(lut,n_tab)
        class(lu_tab_1D) :: lut
        integer, intent(in) :: n_tab
        allocate(lut%table(n_tab))
    end subroutine

    subroutine deallocate_LUT(lut)
        class(lu_tab_1D) :: lut
        deallocate(lut%table)
    end subroutine

    function find_LUT(lut,xtab,x)result(y)
        double precision :: x
        class(lu_tab_1D) :: lut
        type(lu_tab_1D) :: xtab
        double precision :: y

        integer :: n, i
        double precision :: dx, i_float, alpha, eps

        n = size(xtab%table)

        if (n < 2) then
            error stop "Lookup table must have at least two entries"
        end if

        dx  = xtab%table(2) - xtab%table(1)
        eps = abs(dx)*1.0d-12

        ! bounds check (ascending or descending)
        if (x < min(xtab%table(1), xtab%table(n)) - eps .or. &
            x > max(xtab%table(1), xtab%table(n)) + eps) then
            error stop "Requested x is outside lookup table range"
        end if

        ! floating-point index
        i_float = 1.0d0 + (x - xtab%table(1)) / dx

        ! --- exact endpoints ---
        if (i_float <= 1.0d0 + eps) then
            y = lut%table(1)
            return
        end if

        if (i_float >= dble(n) - eps) then
            y = lut%table(n)
            return
        end if

        ! --- interior: interpolate ---
        i = floor(i_float)        ! guaranteed: 1 <= i <= n-1
        alpha = i_float - dble(i)

        y = (1.0d0 - alpha)*lut%table(i) + alpha*lut%table(i+1)

    end function

end module MOD_look_up_table
