! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_work_arrays) SMOD_initialize_work_arrays
    implicit none

contains
    module subroutine initialize_work(work,decomposition)
        use MOD_user_data
        use MOD_user_subroutines
        use MOD_mathematical_functions
        use MOD_work_arrays,only:work_arrays
        use mpidcl,only:decomp_info
        use mpi_f08
        class(work_arrays) :: work
        type(decomp_info), intent(in) :: decomposition

        integer :: n, n_max, i
        double precision :: d_beta
        double precision :: beta

        integer :: ierr,rank

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

        ! Mandatory
        call work%flash_work_to_zero(decomposition)

    end subroutine
end submodule
