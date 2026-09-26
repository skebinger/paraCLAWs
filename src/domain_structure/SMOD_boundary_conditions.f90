! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_domain) SMOD_boundary_conditions
    implicit none

contains

    module subroutine update_rank_boundaries(decomposition,solution)
        !! Performs an update of all rank boundaries.
        !!
        !! First performs halo exchange for all internal boundaries and then applies
        !! physical boundary conditions.
        !! Uses the `mpidcl` library for domain decomposition.
        use mpi_f08
        use mpidcl,only:exchange_halos,get_cartesian_comm,get_neighbouring_ranks,decomp_info
        use MOD_solver_parameters,only:num_equations,num_ghost,num_aux

        type(decomp_info), intent(in) :: decomposition
        type(solution_fields), intent(inout) :: solution

        integer :: size, rank, ierr
        type(MPI_Comm) :: comm_cart
        integer :: nbr_top, nbr_bottom, nbr_left, nbr_right

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)


        ! trigger halo exchange
        ! First: get the cartesian communicator
        comm_cart = decomposition%get_cartesian_comm()
        ! Second retrieve info about neighbouring ranks
        call decomposition%get_neighbouring_ranks(nbr_left, nbr_right, nbr_bottom, nbr_top)

        ! perform the halo exchange
        ! = handle internal boundary conditions; SO FAR I IGNORE PERIODIC BOUNDARIES!!
        call exchange_halos(decomposition,solution%vars,num_equations,num_ghost,nbr_left,nbr_right,nbr_bottom,nbr_top,comm_cart)
        call exchange_halos(decomposition,solution%aux,num_aux,num_ghost,nbr_left,nbr_right,nbr_bottom,nbr_top,comm_cart)
        
        !safe to proceed now without pause since all communications are finished and
        !ranks now operate independently for most of the timestep update
        
        !treat the exterior boundaries
        call apply_physical_boundaries(decomposition,nbr_left,nbr_right,nbr_bottom,nbr_top,solution)

    end subroutine

    subroutine apply_physical_boundaries(decomposition, left, right, bottom, top, solution)
        !! Apply physical boundary treatment on processes touching domain edges
        !!
        !! This subroutine checks if the current process has external (physical) boundaries
        !! and applies boundary conditions on ghost cell layers. Actual BC logic is delegated.
        use mpi_f08,only:MPI_PROC_NULL
        use MOD_solver_parameters,only:num_equations,num_ghost
        use mpidcl,only:decomp_info
        use MOD_user_data
        type(decomp_info), intent(in) :: decomposition
        integer, intent(in) :: left, right, bottom, top !! neighbouring ranks
        type(solution_fields), intent(inout) :: solution

        integer :: i, j, v
        integer :: xi_start, xi_end, eta_start, eta_end

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        ! Left boundary (physical)
        if (left == MPI_PROC_NULL) then
            do j = eta_start, eta_end
                do i = xi_start - num_ghost, xi_start-1
                    error stop "Left boundary condition not implemented yet!"
                end do
            end do
        end if

        ! Right boundary (physical)
        if (right == MPI_PROC_NULL) then
            do j = eta_start, eta_end
                do i = xi_end + 1, xi_end + num_ghost
                    do v = 1, num_equations
                        error stop "Right boundary condition not implemented yet!"
                    end do
                end do
            end do
        end if

        ! Bottom boundary (physical)
        if (bottom == MPI_PROC_NULL) then
            do j = eta_start - num_ghost, eta_start-1
                do i = xi_start, xi_end
                    error stop "Bottom boundary condition not implemented yet!"
                end do
            end do
        end if

        ! Top boundary (physical)
        if (top == MPI_PROC_NULL) then
            do j = eta_end + 1, eta_end + num_ghost
                do i = xi_start, xi_end
                    error stop "Top boundary condition not implemented yet!"
                end do
            end do
        end if

    end subroutine apply_physical_boundaries

end submodule
