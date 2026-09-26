! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

program paraCLAWs
    !! Main program section of the paraCLAWs.

    !parallel libraries:
    use mpi_f08
    use mpidcl,only:decomp_info,initialize_decomposition
    use omp_lib

    !project modules
    use MOD_user_subroutines,only:setup_user_input,set_initial_conditions,setup_auxiliary_quantities,after_step
    use MOD_time_control,only:time_control,advance_snapshot_index,updateTime,get_CFL,get_dt,get_previous_dt
    use MOD_domain,only:update_rank_boundaries,solution_fields,mesh_fields,write_base_info,rescale_coordinates
    use MOD_solve,only:solve_dimensional_split,integrate_source_term,execute_post
    use MOD_utils,only:tic,toc
    use MOD_user_data,only:write_user_data
    use MOD_IO
    use MOD_work_arrays,only:work_arrays,allocate_work_arrays,deallocate_work_arrays

    implicit none
    integer :: exitstat, cmdstat, filestat

    !MPI VARIABLES:
    integer :: ierr, rank, size
    integer :: provided

    !time control info
    type(time_control) :: t_ctrl

    !Decomposition info
    type(decomp_info) :: decomposition

    !Domain arrays
    type(mesh_fields) :: mesh !! container for all mesh fields
    type(solution_fields) :: solution !! container of all solution fields at the current time step
    type(solution_fields) :: solution_old !! container of all solution fields at the previous time step

    !Work arrays
    type(work_arrays) :: work !! container for work arrays

    !===========================================================
    ! Initialize MPI environment
    !===========================================================
    !Initialize MPI as thread safe (only master thread is allowd to make MPI clls)
    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    !define output formats
1   format('paraCLAWs - time: ',f10.8,' current dt: ',f10.8,' CFL(max)=',f10.5,' next dt: ',f10.8,' => snapshot: ', i4)
2   format('paraCLAWs - time: ',f10.8,' current dt: ',f10.8,' CFL(max)=',f10.5,' next dt: ',f10.8)

    if (rank==0)then
        write(*,*) "========================================================="
        call execute_command_line('figlet "paraCLAWs"',.TRUE.,exitstat,cmdstat)
        write(*,*) "parallelised hyperbolic Conservation LaW solver"
        write(*,*) "========================================================="

        call execute_command_line('rm -r output',.TRUE.,exitstat,cmdstat)
        call execute_command_line('mkdir -p output',.TRUE.,exitstat,cmdstat)

        print *, "Number of processes started: ", size
#ifdef _OPENMP
        !$OMP PARALLEL
        !$OMP SINGLE
        print *, "Number of threads per rank: ",omp_get_num_threads()
        !$OMP END SINGLE
        !$OMP END PARALLEL
        write(*,*) "========================================================="
#endif
    end if
    ! Just a small check to see if all ranks started properly
    call MPI_Barrier(MPI_COMM_WORLD)
    write(*,*) "RANK, ", rank, " is ready"
    call MPI_Barrier(MPI_COMM_WORLD)

    !===========================================================
    ! Load solver setup
    !===========================================================
    call read_solver_control(t_ctrl,mesh,solution)

    !===========================================================
    ! Execute the MPI domain decomposition from mpidcl library
    !===========================================================
    call initialize_decomposition(decomposition,mesh%m_xi,mesh%m_eta,0,MPI_COMM_WORLD)
    call decomposition%print_decomposition_summary(MPI_COMM_WORLD)
    call decomposition%setup_cartesian_topology(MPI_COMM_WORLD)
    call decomposition%print_cartesian_rank_layout()

    !===========================================================
    ! Allocate domain memory
    !===========================================================
    call solution%allocate_solution(decomposition)
    call mesh%computational_space%allocate_computational_space(decomposition)
    call mesh%physical_space%allocate_physical_space(decomposition)
    call mesh%quadrilateral_mapping%allocate_quadrilateral_mapping(decomposition)

    !===========================================================
    ! Initialize the computational coordinate system
    !===========================================================
    call mesh%computational_space%initialize_computational_coordinates(decomposition,mesh%m_xi,mesh%m_eta,mesh%xi_dimensions,mesh%eta_dimensions)

    !===========================================================
    ! Load custom user setup
    !===========================================================
    call setup_user_input()
    call rescale_coordinates(mesh,decomposition)

    !===========================================================
    ! Allocate worker arrays
    !===========================================================
    call work%allocate_work_arrays(decomposition)
    ! extra step: might want to do some additional work on the contents of the work object
    call work%initialize_work(decomposition)

    !===========================================================
    ! Setup mesh
    !===========================================================
    call mesh%finalize_mesh(decomposition)

    !===========================================================
    ! Initialize solution field and the aux array
    !===========================================================
    call setup_auxiliary_quantities(decomposition,mesh,solution)
    call set_initial_conditions(decomposition,mesh,solution)

    !force the boundary conditions on the fields to fill the ghost cells
    !at the start time
    call update_rank_boundaries(decomposition,solution)

    !===========================================================
    ! Calculate post process quantities at the initial time
    !===========================================================
    call execute_post(decomposition,mesh,solution)

    !write initial time control data to console
    if(rank==0) write(*,1)t_ctrl%time,t_ctrl%get_previous_dt(),t_ctrl%get_CFL(),t_ctrl%get_dt(),0

    !===========================================================
    ! Write domain settings to disk
    !===========================================================
    call write_base_info(mesh)
    !write mesh to disk
    call write_mesh_to_disk(decomposition,mesh)
    !write initial data to disk
    call writeData(decomposition,t_ctrl%time,solution)
    !repeat writing as a vtk file
    call write_vtk_rank(decomposition,mesh,solution,t_ctrl%time)
    !write decomposition info to disk
    call decomposition%write_decom_to_disk(MPI_COMM_WORLD)
    !write user data to disk
    call write_user_data

    ! Execute the boundary conditions in advance to prepare the
    ! solver for input at t=t_init
    call update_rank_boundaries(decomposition,solution)
    !===========================================================
    ! PERFORM THE MAIN LOOP
    !===========================================================
    do while(t_ctrl%time<t_ctrl%time_final)

        !timestep timer
        if(rank==0) call tic()

        !===========================================================
        ! CREATE A COPY OF THE SOLUTION PRIOR TO ADVANCING
        !===========================================================
        solution_old = solution

        !===========================================================
        ! APPLY THE FRACTIONAL STEP METHOD
        !===========================================================
        !===========================================================
        ! STEP 1: DO THE FLUX DIFFERENCING TO UPDATE Q BASED
        !   ON THE DIMENSIONAL SPLIT METHOD
        !===========================================================
        call solve_dimensional_split(decomposition,t_ctrl,work,mesh,solution,ierr)

        ! Check if the time step was valid
        if(ierr<0)then
            !time step rejected
            if(rank==0) print *, "***************************************************"
            if(rank==0) print *, "WARNING: time step rejected!"

            call t_ctrl%updateTime(step_valid=.false.)

            ! reset the solution
            solution = solution_old

            ! reset the CFL number for the next attempt
            call t_ctrl%reset_CFL()
            ! force the ranks to synchronize
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            !timestep timer
            if(rank==0) call toc()
            if(rank==0) print *, "***************************************************"
            ! skip to the start of the loop
            cycle
        endif

        !===========================================================
        ! STEP 2: INTEGRATE THE SOURCE TERMS ON RHS
        !===========================================================
        call integrate_source_term(decomposition,t_ctrl,mesh,solution)

        !===========================================================
        ! STEP 3: enforce the physical boundary conditions and
        ! update the internal boundaries (halo exchange)
        !===========================================================
        ! Do this here, so that the fields written for post-processing
        ! contain the correct domain boundary values
        call update_rank_boundaries(decomposition,solution)

        !===========================================================
        ! STEP 4: ADVANCE THE TIME BY DT
        !===========================================================
        call t_ctrl%updateTime(step_valid=.true.)

        !===========================================================
        ! HANDLING OF TERMINAL OUTPUT AND WRITING TO DISK
        !===========================================================
        if(t_ctrl%time>=t_ctrl%snapshot_intervals(t_ctrl%get_snapshot_index()))then
            call execute_post(decomposition,mesh,solution)
            !if the current time coincides with an output time defined in time_snapshots
            !=> write to disk
            call writeData(decomposition,t_ctrl%time,solution)
            call write_vtk_rank(decomposition,mesh,solution,t_ctrl%time)
            !call write_postProcessing(current_time_intervals_element)

            if(rank==0) write(*,1) t_ctrl%time,t_ctrl%get_previous_dt(),t_ctrl%get_CFL(),&
                t_ctrl%get_dt(),t_ctrl%get_snapshot_index()
            !increase the counter for the time_intervals by 1
            call t_ctrl%advance_snapshot_index
        else
            !in every other case: just write the current step info to console
            if(rank==0) write(*,2) t_ctrl%time,t_ctrl%get_previous_dt(),t_ctrl%get_CFL(),t_ctrl%get_dt()
        end if

        ! clear the CFL number for the next timestep
        call t_ctrl%reset_CFL()

        call after_step(decomposition,t_ctrl,mesh,solution)

        !Force waiting of all ranks
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !timestep timer
        if(rank==0) call toc()
    end do

    !=========================================================!
    !FINALIZE
    !=========================================================!
    call work%deallocate_work_arrays()
    call mesh%computational_space%deallocate_computational_space()
    call mesh%physical_space%deallocate_physical_space()
    call mesh%quadrilateral_mapping%deallocate_quadrilateral_mapping()
    call solution%deallocate_solution()
    call solution_old%deallocate_solution()
    if(rank==0)then
        print *, "***************************************************"
        print *, "PROGRAM FINISHED SUCCESSFULLY"
        print *, "***************************************************"
    end if

    call MPI_Finalize(ierr)
end program paraCLAWs
