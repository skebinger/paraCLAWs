! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_time_control
    !! This module defines the `time_control` type, which manages simulation time,
    !! timestepping control (including adaptive methods), CFL tracking, and snapshot/output timing.

    implicit none

    !private
    !public :: time_control

    type :: time_control
        !! Object to manage all simulation time tracking, timestep adjustment,
        !! output/snapshot scheduling, and CFL constraints.

        !! Time tracking
        double precision :: time                    !! Current simulation time
        double precision :: time_initial            !! Initial time
        double precision :: time_final              !! Final time

        !! Timestep control
        double precision, private :: dt                !! Current timestep size
        double precision, private :: dt_previous       !! Timestep used in previous iteration
        double precision, private :: dt_max            !! Maximum allowed timestep
        double precision, private :: dt_min            !! Minimum allowed timestep

        !! Snapshot/output timing
        double precision, allocatable :: snapshot_intervals(:) !! Array of times at which output is written to disk
        integer, private :: current_snapshot_index          !! Index of the next snapshot time
        integer, private :: number_of_snapshots             !! Total number of snapshots

        !! CFL (Courant–Friedrichs–Lewy) control
        double precision :: time_precision          !! Precision threshold for timestep adjustment
        double precision, private :: CFL                     !! Maximum CFL encountered this step
        double precision :: CFL_target              !! Target CFL value for adaptive stepping
        double precision :: CFL_previous_timestep   !! CFL value from the previous timestep
        !double precision :: CFL_max_current_step    !! Max CFL recorded in current step
        double precision :: timestep_shrink_limit   !! Minimum factor to shrink timestep
        double precision :: timestep_growth_limit   !! Maximum factor to grow timestep
        double precision, private :: emergency_shrink_ratio = 0.25 !! Emergency ratio to decrease the time step in case of CFL spikes (control rejects time step)

        !! Control flags
        logical :: reset_timestep                   !! Flag to request timestep reset; CURRENTLY UNUSED
        logical :: adaptive_timestepping            !! Enables adaptive timestep control
        logical :: snapshot_exception = .false.     !! Flag to indicate if the current step requires an exception in the adaptive timestep adjustment

    contains
        procedure :: initialize
        procedure :: updateTime
        procedure :: update_max_CFL
        procedure :: advance_snapshot_index
        procedure :: check_timestep_validity

        procedure :: reset_CFL

        ! getter functions:
        procedure :: get_CFL
        procedure :: get_dt
        procedure :: get_previous_dt

        procedure :: get_snapshot_index
    end type time_control

contains

    subroutine initialize(this, time_initial_, time_final_, timestep_, n_output_times_, CFL_target_, &
        timestep_max_, timestep_min_, growth_limit_, shrink_limit_, adaptive_timestepping_)
        !! Initializes the `time_control` object with user-defined simulation and timestep settings.
        !! Also sets up snapshot scheduling and CFL control parameters.

        use mpi_f08
        use MOD_error_handling

        class(time_control), intent(inout) :: this
        double precision, intent(in) :: time_initial_      !! Simulation start time
        double precision, intent(in) :: time_final_        !! Simulation end time
        double precision, intent(in) :: timestep_          !! Initial timestep
        double precision, intent(in) :: timestep_max_      !! Maximum timestep allowed
        double precision, intent(in) :: timestep_min_      !! Minimum timestep allowed
        double precision, intent(in) :: growth_limit_      !! Max growth factor for timestep
        double precision, intent(in) :: shrink_limit_      !! Min shrink factor for timestep
        double precision, intent(in) :: CFL_target_        !! Desired CFL value
        integer, intent(in)          :: n_output_times_    !! Number of snapshot/output times
        logical, intent(in), optional :: adaptive_timestepping_ !! Whether to use adaptive stepping

        double precision :: snapshot_interval
        integer :: i

        integer :: ierr, rank, size

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

        ! Assign time bounds and timestep settings
        this%time_initial = time_initial_
        this%time_final   = time_final_
        this%time         = time_initial_
        this%dt     = timestep_
        this%dt_previous = timestep_
        this%dt_max      = timestep_max_
        this%dt_min      = timestep_min_

        ! Set adaptive timestep flag
        if (present(adaptive_timestepping_)) then
            this%adaptive_timestepping = adaptive_timestepping_
        else
            this%adaptive_timestepping = .False.
        end if

        ! Validate timestep limits
        if (this%dt_min > this%dt_max) then
            if(rank==0) write(*,*) "=============================================================================="
            if(rank==0) write(*,*) "WARNING: Minimum timestep is greater than maximum timestep!"
            if(rank==0) write(*,*) "=============================================================================="
            call end_program("ABORT: Invalid timestep bounds!")
        end if

        ! Ensure snapshot count is at least 1
        if (n_output_times_ < 1) then
            this%number_of_snapshots = 1
            if(rank==0) write(*,*) "=============================================================================="
            if(rank==0) write(*,*) "WARNING: number of output timestamps < 1!"
            if(rank==0) write(*,*) "=> Overwritten; changed to n_output_times = 1"
            if(rank==0) write(*,*) "=============================================================================="
        else
            this%number_of_snapshots = n_output_times_
        end if

        ! Allocate snapshot intervals and compute times
        allocate(this%snapshot_intervals(this%number_of_snapshots))
        snapshot_interval = (this%time_final - this%time_initial) / this%number_of_snapshots

        do i = 1, this%number_of_snapshots
            this%snapshot_intervals(i) = this%time_initial + i * snapshot_interval
        end do
        this%snapshot_intervals(this%number_of_snapshots) = this%time_final
        this%current_snapshot_index = 1

        ! Initialize CFL tracking
        this%CFL                   = 0.d0
        this%CFL_previous_timestep = 0.d0
        this%CFL_target            = CFL_target_
        this%time_precision        = 0.01d0

        ! Store timestep change bounds
        this%timestep_growth_limit = growth_limit_
        this%timestep_shrink_limit = shrink_limit_
    end subroutine initialize

    subroutine updateTime(this,step_valid)
        !! Advances simulation time based on current dt and CFL.
        !! Handles adaptive time stepping and snapshot-aligned time correction.
        use mpi_f08
        use MOD_error_handling,only:end_program

        class(time_control), intent(inout) :: this
        logical, intent(in) :: step_valid

        real(kind=8) :: change_factor
        real(kind=8) :: dt_step !! a temporary copy of the time step width

        integer :: ierr, rank, size

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

        this%dt_previous = this%dt
        dt_step = this%dt

        if(step_valid)then
            ! time step is accepted

            !=======================================================
            ! Adaptive timestep adjustment (if enabled)
            !=======================================================
            if (this%adaptive_timestepping) then
                ! Calculate factor to adjust dt based on CFL deviation from target
                if (this%CFL > 0.d0) then
                    change_factor = this%CFL_target / this%CFL
                else
                    change_factor = 1.d0
                end if

                ! Clamp change factor within user-defined growth ratio
                if (change_factor >= this%timestep_growth_limit) then
                    change_factor = this%timestep_growth_limit
                elseif (abs(change_factor - 1.d0) <= this%time_precision) then
                    change_factor = 1.d0  ! Skip tiny changes
                end if

                ! Apply change
                dt_step = this%dt_previous * change_factor

            end if

            !=======================================================
            ! Align dt with snapshot time if necessary
            !=======================================================
            if (this%time + dt_step > this%snapshot_intervals(this%current_snapshot_index)) then
                dt_step = this%snapshot_intervals(this%current_snapshot_index) - this%time
                this%snapshot_exception = .true.
            else
                this%snapshot_exception = .false.
            end if

        else
            ! The time step is rejected

            ! in case of constant time steps, the program has to be ended, since the
            ! dt can't be changed
            if (.NOT.this%adaptive_timestepping) then
                ! gracefully shut down the program
                call end_program(msg='constant time steps do not allow for adaptive correction to severe CFL violation!')
            end if

            ! Now actually decrease the time step
            dt_step = this%emergency_shrink_ratio*this%dt
            if(rank==0) print *, "BAD CFL number ", this%CFL
            if(rank==0) print *, "-> retry with dt=", dt_step
        end if

        ! Enforce min/max dt
        if (dt_step > this%dt_max) dt_step = this%dt_max
        if (dt_step < this%dt_min) then
            if(rank==0) write(*,*) "WARNING: dt limited to dt_min =", this%dt_min
            if (.NOT.this%adaptive_timestepping.AND.(.NOT.this%snapshot_exception)) then
                ! gracefully shut down the program
                call end_program(msg='minimum time step size reached, but control still detects CFL violation - constant dt not compatible!')
            end if
            dt_step = this%dt_min
        end if

        !=======================================================
        ! Apply the time step
        !=======================================================
        if(step_valid)then
            this%time = this%time + dt_step
        end if
        this%dt=dt_step

        ! Track CFL values
        !this%CFL_previous_timestep = this%CFL
        !this%CFL_max_current_step = this%CFL

        ! Reset CFL accumulator
        !this%CFL = 0.d0
    end subroutine updateTime

    subroutine update_max_CFL(this, CFL_curr)
        !! Updates the maximum CFL number tracked for the current timestep.

        class(time_control), intent(inout) :: this
        double precision, intent(in) :: CFL_curr  !! Current CFL value to compare against

        this%CFL = max(this%CFL, CFL_curr)
    end subroutine update_max_CFL

    subroutine advance_snapshot_index(this)
        !! Advances the index for the next output snapshot time.

        class(time_control), intent(inout) :: this

        this%current_snapshot_index = this%current_snapshot_index + 1
    end subroutine advance_snapshot_index

    function check_timestep_validity(this)result(valid_flag)
        class(time_control) :: this
        logical:: valid_flag !! flag determining if current time-step is valid

        ! Check if the time step would need to shrink by more than the allowed limit
        valid_flag = this%CFL_target/this%CFL>=this%timestep_shrink_limit
    end function

    subroutine reset_CFL(this)
        !! Resets the CFL tracking for the current timestep.
        !! This should be called at the end of each time step after CFL values have been used for timestep adjustment.
        !! Note: The previous timestep's CFL is stored before resetting, allowing for diagnostics or logging.

        class(time_control), intent(inout) :: this

        this%CFL_previous_timestep = this%CFL
        this%CFL = 0.d0
    end subroutine reset_CFL

    function get_CFL(this) result(CFL_out)
        !! Returns the current CFL number for this time step.

        class(time_control) :: this
        double precision :: CFL_out

        CFL_out = this%CFL
    end function get_CFL

    function get_dt(this) result(dt_out)
        !! Returns the current time step size.

        class(time_control) :: this
        double precision :: dt_out

        dt_out = this%dt
    end function get_dt

    function get_previous_dt(this) result(dt_prev_out)
        !! Returns the time step size from the previous iteration.

        class(time_control) :: this
        double precision :: dt_prev_out

        dt_prev_out = this%dt_previous
    end function get_previous_dt

    function get_snapshot_index(this) result(index_out)
        !! Returns the index of the next snapshot time.

        class(time_control) :: this
        integer :: index_out

        index_out = this%current_snapshot_index
    end function get_snapshot_index

end module MOD_time_control
