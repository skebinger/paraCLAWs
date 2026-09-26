! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_user_data
    !! Module to store persistent specific user data

    implicit none

contains

    subroutine write_user_data()
        use mpi_f08

        integer :: size,rank,ierr

        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        if(rank==0)then
            open(unit=4,file='output/user_data.dat',status='unknown',form='formatted')

            close(unit=4)
        end if

    end subroutine

end module
