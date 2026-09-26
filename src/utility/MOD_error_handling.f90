! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_error_handling
    use mpi_f08
    implicit none

    abstract interface
        subroutine error_callback()
        end subroutine
    end interface
contains

    subroutine end_program(msg)
        !! Used for a premature, but graceful shutdown of the program
        character(len=*), intent(in) :: msg
        integer :: rank
        integer :: ierr

        call MPI_Finalize(ierr)

        print *, "***************************************************"
        print *, "PROGRAM FINISHED PREMATURELY"
        print *, "REASON: "
        print *, msg
        print *, "***************************************************"

        stop
    end subroutine

    subroutine abort_program(msg, comm)
        !! Used if a fatal error is encountered by one of the ranks
        !! and a graceful shutdown is too tedious.
        !! For example if too deep into a subroutine hirarchy.
        character(len=*), intent(in) :: msg
        type(MPI_Comm), intent(in) :: comm
        integer :: rank
        integer :: ierr

        call MPI_Comm_rank(comm, rank, ierr)

        write(*,*) 'Rank ', rank, ' encountered error: ', trim(msg)
        call MPI_Abort(comm, 1, ierr)
    end subroutine abort_program

end module MOD_error_handling
