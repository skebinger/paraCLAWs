! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_solve
    !! Module containing subroutines to advance a set of hyperbolic PDEs in time.
    !!
    !! Mainly the containing methods form an algorithm performing the fractional
    !! step method. Currently, the sweeps are structured with a dimensional split method
    !! and the source (right hand side) is evaluated last.
    implicit none

    private :: integrate_flux_difference
    public :: solve_dimensional_split

contains

    subroutine solve_dimensional_split(decomposition,t_ctrl,work,mesh,solution,ierr)
        !! Subroutine calculating the homogeneous solution for a set of
        !! hyperbolic PDEs with the dimensional split method.
        !!
        !! Decomposing the jump at i-1/2 or j-1/2 and the integration are in one loop but
        !! executed out of sync. To reduce nested loops, the integration is delayed to allow
        !! calculation of the waves and eigenvalues at the interface and then executed after the
        !! delay is passed.

        !! error return codes:
        !!  1 = success
        !!  -1 = time step rejected
        use MOD_time_control,only:time_control
        use mpidcl,only:decomp_info,find_min_max_scalar
        use MOD_solver_parameters,only:num_waves,add_higher_order_flux_correction
        use MOD_flux_calculation,only:wave_decomposition,split_flux
        use MOD_domain,only:solution_fields,mesh_fields
        use MOD_flux_correction,only:limiter,second_order_flux_correction
        use MOD_work_arrays,only:work_arrays,flash_work_to_zero
        use mpi_f08

        type(decomp_info), intent(in) :: decomposition
        type(time_control), intent(inout) :: t_ctrl
        type(work_arrays), intent(inout) :: work !! container for work arrays
        type(mesh_fields) :: mesh !! container for all mesh fields
        type(solution_fields) :: solution !! container of all solution fields at the current time step
        integer, intent(inout) :: ierr !! an error code for signalling to the calling subroutine

        integer :: ixy ! current sweep direction
        integer :: xi_start, xi_end, eta_start, eta_end
        integer :: i,j,mw
        integer :: local_validity, global_validity
        logical :: step_valid

        double precision :: CFL_local
        double precision :: CFL_max_all

        CFL_local = 0
        CFL_max_all = 0

        !! assign the error code (just assume failure for now)
        ierr = -1

        ! get MPI block bounds
        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        !===================================================================
        ! xi-sweep
        ixy=1
        !$OMP PARALLEL
        !$OMP DO PRIVATE(i,mw) REDUCTION(MAX:CFL_local)
        do j=eta_start,eta_end !slice direction
            do i=xi_start-1,xi_end+1 !integration direction
                ! Calculate eigenvalues and waves
                call wave_decomposition(ixy,i,j,mesh,solution,work,decomposition)

                ! update max CFL
                do mw=1,num_waves
                    CFL_local = max(CFL_local,work%lambda(mw,i,j)*t_ctrl%get_dt()/mesh%computational_space%d_xi)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        call split_flux(decomposition,work%lambda,work%waves,work%Aminus_dQ,work%Aplus_dQ)

        call integrate_flux_difference(decomposition,t_ctrl,work,mesh,solution,ixy)

        ! these two nested loops may be combined with the first sweep! just add another delay to loop counter
        if(add_higher_order_flux_correction.eqv..true.)then
            call limiter(decomposition,work,ixy)

            call second_order_flux_correction(decomposition,t_ctrl,work,mesh,ixy)

            call integrator_flux_correction_increment(decomposition,t_ctrl,work,mesh,solution,ixy)
        end if

        ! identify and communicate the (possibly new) max CFL to all ranks
        call find_min_max_scalar(CFL_local,CFL_max_all,MPI_MAX,MPI_COMM_WORLD)
        call t_ctrl%update_max_CFL(CFL_max_all)
        !===================================================================
        call work%flash_work_to_zero(decomposition)

        !===================================================================
        ! check if sweep is valid (CFL criterion not violated)
        step_valid = t_ctrl%check_timestep_validity()
        local_validity = merge(1, 0, step_valid)

        !===================================================================
        ! SYNCHRONIZE VALIDITY DECISION ACROSS ALL RANKS AFTER XI-SWEEP
        ! If ANY rank rejects after the first sweep, all ranks can abort before wasting
        ! time on the eta-sweep.
        call find_min_max_scalar(local_validity,global_validity,MPI_MIN,MPI_COMM_WORLD)
        if(global_validity == 0) then
            ierr = -1
            return
        endif

        !===================================================================
        ! Continue to eta-sweep
        ixy=2
        ! possible to switch i and j indizes?
        !$OMP PARALLEL
        !$OMP DO PRIVATE(j,mw) REDUCTION(MAX:CFL_local)
        do i=xi_start,xi_end
            do j=eta_start-1,eta_end+1
                ! Calculate eigenvalues and waves
                call wave_decomposition(ixy,i,j,mesh,solution,work,decomposition)
                ! update max CFL
                do mw=1,num_waves
                    CFL_local=max(CFL_local,work%lambda(mw,i,j)*t_ctrl%get_dt()/mesh%computational_space%d_eta)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        call split_flux(decomposition,work%lambda,work%waves,work%Aminus_dQ,work%Aplus_dQ)

        call integrate_flux_difference(decomposition,t_ctrl,work,mesh,solution,ixy)

        if(add_higher_order_flux_correction.eqv..true.)then
            call limiter(decomposition,work,ixy)

            call second_order_flux_correction(decomposition,t_ctrl,work,mesh,ixy)

            call integrator_flux_correction_increment(decomposition,t_ctrl,work,mesh,solution,ixy)
        end if

        ! identify and communicate the (possibly new) max CFL to all ranks
        call find_min_max_scalar(CFL_local,CFL_max_all,MPI_MAX,MPI_COMM_WORLD)
        call t_ctrl%update_max_CFL(CFL_max_all)
        !===================================================================
        call work%flash_work_to_zero(decomposition)

        !===================================================================
        ! check if sweep is valid (CFL criterion not violated)
        step_valid = t_ctrl%check_timestep_validity()
        ! Keep local_validity as 0 if xi-sweep failed (using AND logic)
        local_validity = local_validity * merge(1, 0, step_valid)

        !===================================================================
        ! SYNCHRONIZE VALIDITY DECISION ACROSS ALL RANKS
        ! If ANY rank rejects, ALL ranks must reject to avoid MPI deadlock
        call find_min_max_scalar(local_validity,global_validity,MPI_MIN,MPI_COMM_WORLD)

        if(global_validity == 0)then
            ierr = -1
        else
            ierr = 1
        endif
    end subroutine

    subroutine integrate_flux_difference(decomposition,t_ctrl,work,mesh,solution,ixy)
        !! Integrates the flux difference for each component $q_p$ at the cell i,j
        !! $$ q_{p,i,j}^{n+1} = q_{p,i,j}^n - \frac{\Delta t}{\Delta x \kappa_{i,j}} (\gamma_{i-1/2} A^+\Delta Q _{p,i-1/2} + \gamma_{i+1/2} A^-\Delta Q_{p,i+1/2}) $$
        !! Considering the grid mapping, the fluxes are rescaled with the length ratios $\gamma_{i\pm 1/2}$ to ensure the correct physical fluxes are integrated. <br>
        !! The calculation for the eta-sweep is analogous.
        use MOD_solver_parameters,only:num_equations
        use mpidcl,only:decomp_info
        use MOD_time_control,only:time_control
        use MOD_work_arrays,only:work_arrays
        use MOD_domain,only:mesh_fields,solution_fields
        type(decomp_info), intent(in) :: decomposition
        type(time_control), intent(in) :: t_ctrl
        type(work_arrays), intent(in) :: work
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution
        integer, intent(in) :: ixy

        integer :: me
        integer :: i,j
        integer :: xi_start, xi_end, eta_start, eta_end
        double precision :: delta, dt_over_delta
        double precision :: len_rat_l,len_rat_r

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        delta = merge(mesh%computational_space%d_xi, mesh%computational_space%d_eta, ixy==1)

        dt_over_delta = t_ctrl%get_dt()/delta

        ! cell centred update
        if(ixy==1)then
            !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,me)
            do j=eta_start,eta_end
                do i=xi_start,xi_end
                    ! determine the rescaling of the fluxes due to grid mapping
                    len_rat_l = mesh%quadrilateral_mapping%length_ratio_xi(1,i,j)
                    len_rat_r = mesh%quadrilateral_mapping%length_ratio_xi(1,i+1,j)
                    do me=1,num_equations
                        solution%vars(me,i,j) = solution%vars(me,i,j) - &
                            dt_over_delta/mesh%quadrilateral_mapping%capacity(i,j) * &
                            (len_rat_l*work%Aplus_dQ (me,i,j) + len_rat_r*work%Aminus_dQ (me,i+1,j))
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,me)
            do j=eta_start,eta_end
                do i=xi_start,xi_end
                    ! determine the rescaling of the fluxes due to grid mapping
                    len_rat_l = mesh%quadrilateral_mapping%length_ratio_eta(1,i,j)
                    len_rat_r = mesh%quadrilateral_mapping%length_ratio_eta(1,i,j+1)
                    do me=1,num_equations
                        solution%vars(me,i,j) = solution%vars(me,i,j) - &
                            dt_over_delta/mesh%quadrilateral_mapping%capacity(i,j) * &
                            (len_rat_l*work%Aplus_dQ (me,i,j) + len_rat_r*work%Aminus_dQ (me,i,j+1))
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if
    end subroutine

    subroutine integrator_flux_correction_increment(decomposition,t_ctrl,work,mesh,solution,ixy)
        !! Subroutine to add the flux correction increment to the solution after the second order flux correction is calculated.
        !! $$ q_{p,i,j}^{n+1} = q_{p,i,j}^{n+1} - \frac{\Delta t}{\Delta x \kappa_{i,j}} (f_{p,i+1,j}-f_{p,i,j}) $$
        !! The calculation for the eta-sweep is analogous.
        use MOD_solver_parameters,only:num_equations
        use mpidcl,only:decomp_info
        use MOD_time_control,only:time_control
        use MOD_work_arrays,only:work_arrays
        use MOD_domain,only:mesh_fields,solution_fields
        type(decomp_info), intent(in) :: decomposition
        type(time_control), intent(in) :: t_ctrl
        type(work_arrays), intent(in) :: work
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution
        integer, intent(in) :: ixy

        integer :: me
        integer :: i,j
        integer :: xi_start, xi_end, eta_start, eta_end
        double precision :: delta, dt_over_delta

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        delta = merge(mesh%computational_space%d_xi, mesh%computational_space%d_eta, ixy==1)

        dt_over_delta = t_ctrl%get_dt()/delta

        ! cell centred update
        if(ixy==1)then
            !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,me)
            do j=eta_start,eta_end
                do i=xi_start,xi_end
                    do me=1,num_equations
                        solution%vars(me,i,j) = solution%vars(me,i,j) - &
                            dt_over_delta/mesh%quadrilateral_mapping%capacity(i,j) * &
                            (work%f(me,i+1,j)-work%f(me,i,j))
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,me)
            do j=eta_start,eta_end
                do i=xi_start,xi_end
                    do me=1,num_equations
                        solution%vars(me,i,j) = solution%vars(me,i,j) - &
                            dt_over_delta/mesh%quadrilateral_mapping%capacity(i,j) * &
                            (work%f(me,i,j+1)-work%f(me,i,j))
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if
    end subroutine

    subroutine integrate_source_term(decomposition,t_ctrl,mesh,solution)
        !! Integrates the source term for each component $q_p$ at the cell i,j according to the fractional step method.
        !! $$ q_{p,i,j}^{n+1} = q_{p,i,j}^n + \phi_p \Delta t $$
        !! Currently implements a simple (forward) Euler method for the source term integration.
        use MOD_solver_parameters,only:integrate_source,num_equations
        use mpidcl,only:decomp_info
        use MOD_time_control,only: time_control
        use MOD_domain,only:mesh_fields,solution_fields
        use MOD_user_subroutines,only:source_function

        type(decomp_info), intent(in) :: decomposition
        type(time_control), intent(in) :: t_ctrl
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution

        integer :: me
        integer :: xi_start, xi_end, eta_start, eta_end
        integer :: i,j
        double precision :: phi(num_equations)

        if(integrate_source)then

            call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

            !$OMP PARALLEL DO PRIVATE(i,me,phi)
            do j=eta_start,eta_end
                do i=xi_start,xi_end

                    phi = source_function(i,j,mesh,solution,decomposition)

                    do me=1,num_equations

                        solution%vars(me,i,j) = solution%vars(me,i,j) + phi(me)*t_ctrl%get_dt()

                    end do
                end do
            end do
            !$OMP END PARALLEL DO

        end if

    end subroutine

    subroutine execute_post(decomposition,mesh,solution)
        !! Subroutine executing the post processing.
        !! Calls the post subroutine to be filled by user cell by cell.
        use MOD_time_control,only: time_control
        use mpidcl,only:decomp_info
        use MOD_user_subroutines,only:user_post
        use MOD_domain

        type(decomp_info), intent(in) :: decomposition
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution

        integer :: xi_start, xi_end, eta_start, eta_end
        integer :: i,j

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        !$OMP PARALLEL DO PRIVATE(i)
        do j=eta_start,eta_end
            do i=xi_start,xi_end
                call user_post(i,j,mesh,solution,decomposition)
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine

end module MOD_solve
