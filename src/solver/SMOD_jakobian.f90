! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule (MOD_flux_calculation) SMOD_jakobian
    !!Submodule containing jakobian related calculations
    implicit none

contains

    module function getJakobi(ixy,i,j,mesh,solution,work,decomposition)result(jakobi)
        !! Returns the jakobian matrix for the state vars(i,j); w/ aux array
        use mpi_f08
        use MOD_solver_parameters,only:num_equations,psi=>grav
        use MOD_user_subroutines
        use MOD_user_data
        use MOD_index
        use MOD_domain,only:solution_fields,mesh_fields
        use MOD_mathematical_functions, only: integrate
        use MOD_work_arrays,only:work_arrays
        use mpidcl,only:decomp_info
        use MOD_error_handling,only:abort_program
        integer, intent(in) :: ixy !! direction of sweep
        integer, intent(in) :: i,j
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(in) :: solution
        type(work_arrays), intent(in) :: work
        type(decomp_info), intent(in) :: decomposition
        double precision :: jakobi(num_equations,num_equations)

        integer :: k,l
        integer :: rank,ierr
        !==========================================================
        ! START OF USER INPUT
        !==========================================================

        call abort_program(msg="getJakobi is not implemented yet. Please implement it in the submodule SMOD_jakobian",comm=MPI_COMM_WORLD)

        do k=1,num_equations
            do l=1,num_equations
                if(isnan(jakobi(k,l)))then
                    write(*,*)"isnan in jakobi"
                    write(*,*) "cell: ", i,j
                    write(*,*) "direction: ", ixy
                    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
                    print *, "rank ", rank
                    write(*,*) solution%vars(:,i,j)
                    write(*,*) jakobi(1,:),1
                    write(*,*) jakobi(2,:),2
                    write(*,*) jakobi(3,:),3
                    error stop
                end if
            end do
        end do
    end function

end submodule
