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

    integer, parameter :: BC_USER=1, BC_ZERO_GRADIENT=2, BC_CONSTANT_GRADIENT=3

contains

    module subroutine update_rank_boundaries(decomposition,mesh,solution)
        !! Performs an update of all rank boundaries.
        !!
        !! First performs halo exchange for all internal boundaries and then applies
        !! physical boundary conditions.
        !! Uses the `mpidcl` library for domain decomposition.
        use mpi_f08
        use mpidcl,only:exchange_halos,get_cartesian_comm,get_neighbouring_ranks,decomp_info
        use MOD_solver_parameters,only:num_equations,num_ghost,num_aux

        type(decomp_info), intent(in) :: decomposition
        type(mesh_fields), intent(in) :: mesh
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
        call apply_physical_boundaries(decomposition,nbr_left,nbr_right,nbr_bottom,nbr_top,mesh,solution)

    end subroutine

    subroutine apply_physical_boundaries(decomposition, left, right, bottom, top, mesh, solution)
        !! Apply physical boundary treatment on processes touching domain edges
        !!
        !! This subroutine checks if the current process has external (physical) boundaries
        !! and applies boundary conditions on ghost cell layers.
        use mpi_f08,only:MPI_PROC_NULL
        use MOD_solver_parameters,only:num_equations,num_ghost
        use mpidcl,only:decomp_info
        use MOD_user_data
        use MOD_error_handling
        type(decomp_info), intent(in) :: decomposition
        integer, intent(in) :: left, right, bottom, top !! neighbouring ranks
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution

        integer :: i, j, me
        integer :: ilow, ihigh, jlow, jhigh

        integer :: i_reverse ! index to count in reverse direction
        double precision :: width, delta_var

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        ! Left boundary (physical)
        if (left == MPI_PROC_NULL) then
            do j = jlow, jhigh
                do i = ilow - num_ghost, ilow-1

                    select case(mesh%bc_xi_low)
                      case(BC_USER)
                        ! Apply user-defined boundary condition
                        call abort_program(msg="User-defined boundary condition at lower xi boundary not implemented yet!", comm=MPI_COMM_WORLD)

                      case(BC_ZERO_GRADIENT)
                        ! Apply zero gradient boundary condition
                        do me = 1, num_equations
                            solution%vars(me,i,j) = solution%vars(me,ilow,j)
                        end do

                      case(BC_CONSTANT_GRADIENT)
                        ! Linearly extrapolate from the two nearest cells.
                        call abort_program(msg="Constant gradient boundary condition not implemented yet!", comm=MPI_COMM_WORLD)

                      case default
                        call abort_program(msg="Left boundary condition not implemented yet!", comm=MPI_COMM_WORLD)
                    end select

                end do
            end do
        end if

        ! Right boundary (physical)
        if (right == MPI_PROC_NULL) then
            do j = jlow, jhigh
                do i = ihigh + 1, ihigh + num_ghost

                    select case(mesh%bc_xi_low)
                      case(BC_USER)
                        ! Apply user-defined boundary condition
                        call abort_program(msg="User-defined boundary condition at upper xi boundary not implemented yet!", comm=MPI_COMM_WORLD)

                      case(BC_ZERO_GRADIENT)
                        ! Apply zero gradient boundary condition
                        do me = 1, num_equations
                            solution%vars(me,i,j) = solution%vars(me,ihigh,j)
                        end do

                      case(BC_CONSTANT_GRADIENT)
                        ! Linearly extrapolate from the two nearest cells.
                        call abort_program(msg="Constant gradient boundary condition not implemented yet!", comm=MPI_COMM_WORLD)

                      case default
                        call abort_program(msg="Left boundary condition not implemented yet!", comm=MPI_COMM_WORLD)
                    end select
                end do
            end do
        end if

        ! Bottom boundary (physical)
        if (bottom == MPI_PROC_NULL) then
            do j = jlow - num_ghost, jlow-1
                do i = ilow, ihigh

                    select case(mesh%bc_xi_low)
                      case(BC_USER)
                        ! Apply user-defined boundary condition
                        call abort_program(msg="User-defined boundary condition at lower eta boundary not implemented yet!", comm=MPI_COMM_WORLD)

                      case(BC_ZERO_GRADIENT)
                        ! Apply zero gradient boundary condition
                        do me = 1, num_equations
                            solution%vars(me,i,j) = solution%vars(me,i,jlow)
                        end do

                      case(BC_CONSTANT_GRADIENT)
                        ! Linearly extrapolate from the two nearest cells.
                        call abort_program(msg="Constant gradient boundary condition not implemented yet!", comm=MPI_COMM_WORLD)

                      case default
                        call abort_program(msg="Bottom boundary condition not implemented yet!", comm=MPI_COMM_WORLD)
                    end select

                end do
            end do
        end if

        ! Top boundary (physical)
        if (top == MPI_PROC_NULL) then
            do j = jhigh + 1, jhigh + num_ghost
                do i = ilow, ihigh

                    select case(mesh%bc_xi_high)
                      case(BC_USER)
                        ! Apply user-defined boundary condition
                        call abort_program(msg="User-defined boundary condition at upper xi boundary not implemented yet!", comm=MPI_COMM_WORLD)

                      case(BC_ZERO_GRADIENT)
                        ! Apply zero gradient boundary condition
                        do me = 1, num_equations
                            solution%vars(me,i,j) = solution%vars(me,i,jhigh)
                        end do

                      case(BC_CONSTANT_GRADIENT)
                        ! Linearly extrapolate from the two nearest cells.
                        call abort_program(msg="Constant gradient boundary condition not implemented yet!", comm=MPI_COMM_WORLD)

                      case default
                        call abort_program(msg="Top boundary condition not implemented yet!", comm=MPI_COMM_WORLD)
                    end select

                end do
            end do
        end if

    end subroutine apply_physical_boundaries

end submodule
