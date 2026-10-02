! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_user_subroutines) SMOD_set_initial_conditions
    implicit none
contains

    module subroutine set_initial_conditions(decomposition,mesh,solution)
        use mpidcl,only:decomp_info
        use MOD_domain,only:solution_fields,mesh_fields
        use MOD_solver_parameters,only:num_ghost
        use MOD_user_data
        use MOD_index
        use mpi_f08

        type(decomp_info), intent(in) :: decomposition
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution

        integer :: i,j
        integer :: xi_start, xi_end, eta_start, eta_end
        integer :: rank,ierr

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        !$OMP PARALLEL DO PRIVATE(i)
        do j=eta_start-num_ghost,eta_end+num_ghost
            do i=xi_start-num_ghost,xi_end+num_ghost
                if(sqrt(mesh%physical_space%x(i,j)**2 + mesh%physical_space%y(i,j)**2) < 0.5d0)then
                    solution%vars(1,i,j) = 2.d0
                else
                    solution%vars(1,i,j) = 1d0
                end if
                solution%vars(2,i,j) = 0.0d0
                solution%vars(3,i,j) = 0.0d0
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine

end submodule
