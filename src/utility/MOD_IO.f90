! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !


module MOD_IO
    !! Module containing the routines for reading and writing data to disk.
    !!
    !! Every instance of writing, each rank writes its own data into a rank specific directory.
    implicit none

    private :: create_directory
#ifdef USE_VTK
    private :: pvd_times, write_pvd_collection
    double precision, allocatable, save :: pvd_times(:)
#endif

contains

    subroutine read_solver_control(t_ctrl,mesh,solution)
        !! Subroutine to read the solver control parameters from the control.in file.
        use MOD_solver_parameters,only:set_solver_parameters
        use MOD_time_control,only:time_control
        use MOD_domain,only:mesh_fields,solution_fields

        type(time_control), intent(inout) :: t_ctrl
        type(mesh_fields), intent(inout) :: mesh
        type(solution_fields), intent(inout) :: solution

        integer :: filestat
        integer :: num_equations
        integer :: num_dimensions
        integer :: num_ghost
        integer :: num_aux
        integer :: num_waves
        double precision :: t_start, t_final
        integer :: number_of_snapshots
        double precision :: dt_initial
        logical :: adaptive_timestepping
        double precision :: CFL_limit
        double precision :: dt_min, dt_max
        double precision :: dt_max_growth_ratio, dt_min_shrink_ratio
        logical :: higher_order_flux_correction
        logical :: use_limiter
        integer :: limiter_method
        integer :: m_xi, m_eta
        double precision :: xi_dimensions(2), eta_dimensions(2)
        integer :: bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper

        open(unit=1, file='control.in', status='old', action='read', iostat=filestat)

        if (filestat /= 0) then
            error stop "Error opening <control.in>. Please ensure the file exists and is in the correct location."
        else
            ! read the general solver parameters to determine the array dimensions
            read(1,*) num_equations
            read(1,*) num_dimensions
            read(1,*) num_ghost
            read(1,*) num_aux
            read(1,*) num_waves

            ! read the time data from the control file
            read(1,*) t_start
            read(1,*) t_final
            read(1,*) number_of_snapshots
            read(1,*) dt_initial
            read(1,*) adaptive_timestepping
            read(1,*) CFL_limit
            read(1,*) dt_min
            read(1,*) dt_max
            read(1,*) dt_max_growth_ratio
            read(1,*) dt_min_shrink_ratio
            read(1,*) higher_order_flux_correction
            read(1,*) use_limiter
            read(1,*) limiter_method

            ! read the mesh data from the control file
            read(1,*) m_xi
            read(1,*) m_eta
            read(1,*) xi_dimensions(1)
            read(1,*) xi_dimensions(2)
            read(1,*) eta_dimensions(1)
            read(1,*) eta_dimensions(2)

            ! read the boundary conditions
            read(1,*) bc_xi_lower
            read(1,*) bc_xi_upper
            read(1,*) bc_eta_lower
            read(1,*) bc_eta_upper
        endif

        close(1)

        ! set the solver control parameters
        call set_solver_parameters(num_equations,num_waves,num_ghost,num_aux,higher_order_flux_correction,use_limiter,limiter_method)

        ! initialize the time control
        call t_ctrl%initialize(t_start, t_final, dt_initial, number_of_snapshots, CFL_limit, &
            dt_max, dt_min, dt_max_growth_ratio, dt_min_shrink_ratio, adaptive_timestepping)

        ! initialize the base mesh type
        call mesh%setup_mesh_base(m_xi, m_eta, bc_xi_lower, bc_xi_upper, bc_eta_lower, bc_eta_upper, xi_dimensions, eta_dimensions)

    end subroutine read_solver_control

    subroutine write_vtk_rank(info,mesh,solution,time)
        use mpi_f08
        use mpidcl,only:decomp_info
        use MOD_domain
        use MOD_solver_parameters,only:num_equations,num_ghost
#ifdef USE_VTK
        use vtk_fortran, only: vtk_file
        use penf,        only: I4P, R8P
#endif

        implicit none

        type(decomp_info), intent(in) :: info
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(in) :: solution
        double precision, intent(in) :: time

        integer :: rank,ierr
        integer :: ilow,ihigh,jlow,jhigh
        integer :: vtk_ilow,vtk_ihigh,vtk_jlow,vtk_jhigh
        integer :: nbr_left,nbr_right,nbr_bottom,nbr_top
        integer :: nx,ny,npts
        integer :: i,j,k,m

        character(len=200) :: filename
        character(len=40) :: rank_str
        character(len=40) :: dirname,subdirname,timestamp_str
        character(len=40) :: q_name
        integer :: error

        ! Unless VTK output is enabled when compiling, this subroutine does nothing
#ifdef USE_VTK
        type(vtk_file) :: vtkobj
        real(R8P), allocatable :: xpts(:),ypts(:),zpts(:)
        real(R8P), allocatable :: qfield(:)

        ! ----------------------------------------------------------------------
        ! Get MPI rank and local mesh bounds
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
        call info%get_local_block_bounds(ilow,ihigh,jlow,jhigh)
        call info%get_neighbouring_ranks(left=nbr_left,right=nbr_right,bottom=nbr_bottom,top=nbr_top)

        ! Include the lower-index neighbor's final point at the start of each piece.
        vtk_ilow = ilow
        vtk_ihigh = ihigh
        vtk_jlow = jlow
        vtk_jhigh = jhigh
        if(nbr_left/=MPI_PROC_NULL) vtk_ilow = ilow - 1
        if(nbr_bottom/=MPI_PROC_NULL) vtk_jlow = jlow - 1
        if((vtk_ilow<ilow .or. vtk_jlow<jlow) .and. num_ghost<1) then
            write(*,*) "Parallel VTK output requires at least one ghost layer."
            call MPI_Abort(MPI_COMM_WORLD,1,ierr)
        end if

        write(rank_str,'(I0)') rank

        nx = vtk_ihigh - vtk_ilow + 1
        ny = vtk_jhigh - vtk_jlow + 1
        npts = nx * ny

        allocate(xpts(npts), ypts(npts), zpts(npts))

        ! ----------------------------------------------------------------------
        ! Flatten mesh coordinates (2D: z = 0)
        k = 0
        do j = vtk_jlow, vtk_jhigh
            do i = vtk_ilow, vtk_ihigh
                k = k + 1
                xpts(k) = mesh%physical_space%x(i,j)
                ypts(k) = mesh%physical_space%y(i,j)
                zpts(k) = 0.0_R8P
            end do
        end do

        ! ----------------------------------------------------------------------
        ! Setup output directories matching writeData style
        write(timestamp_str, '(F20.10)') time
        dirname = 'output/' // trim(adjustl(timestamp_str))
        call create_directory(dirname, rank)
        subdirname = adjustl(trim(dirname) // '/rank_' // trim(rank_str))
        call create_directory(subdirname, rank)

        ! Setup filename for this rank in current timestep's directory
        filename = trim(subdirname) // '/solution_rank_' // trim(rank_str) // '.vts'

        ! Initialize VTK file
        error = vtkobj%initialize(format='binary', filename=filename, &
            mesh_topology='StructuredGrid', &
            nx1=vtk_ilow, nx2=vtk_ihigh, ny1=vtk_jlow, ny2=vtk_jhigh, nz1=1, nz2=1)
        if (error /= 0) stop "VTK initialize failed"

        ! Write geometry
        error = vtkobj%xml_writer%write_piece(nx1=vtk_ilow, nx2=vtk_ihigh, ny1=vtk_jlow, ny2=vtk_jhigh, nz1=1, nz2=1)
        error = vtkobj%xml_writer%write_geo(n=npts, x=xpts, y=ypts, z=zpts)

        ! ----------------------------------------------------------------------
        ! Write solution fields
        error = vtkobj%xml_writer%write_dataarray(location='node', action='open')
        if (error /= 0) stop "VTK write_dataarray open failed"

        do m = 1, num_equations
            allocate(qfield(npts))
            k = 0
            do j = vtk_jlow, vtk_jhigh
                do i = vtk_ilow, vtk_ihigh
                    k = k + 1
                    qfield(k) = solution%vars(m,i,j)
                end do
            end do
            write(q_name, '(A,I0)') 'q', m
            error = vtkobj%xml_writer%write_dataarray(data_name=trim(q_name), x=qfield)
            if (error /= 0) stop "VTK write_dataarray data failed"
            deallocate(qfield)
        end do

        error = vtkobj%xml_writer%write_dataarray(location='node', action='close')
        if (error /= 0) stop "VTK write_dataarray close failed"

        ! Close piece tag
        error = vtkobj%xml_writer%write_piece()
        if (error /= 0) stop "VTK write_piece close failed"

        ! Finalize the VTK file
        error = vtkobj%finalize()

        deallocate(xpts, ypts, zpts)
#endif
    end subroutine

    subroutine write_pvts_rank(info,time)
        use mpi_f08
        use mpidcl,only:decomp_info
        use MOD_solver_parameters,only:num_equations,num_ghost
#ifdef USE_VTK
        use vtk_fortran,only:pvtk_file
        use penf,only:I4P
#endif

        type(decomp_info), intent(in) :: info
        double precision, intent(in) :: time

#ifdef USE_VTK
        type(pvtk_file) :: pvts
        integer :: rank,size,ierr,error
        integer :: ilow,ihigh,jlow,jhigh
        integer :: vtk_ilow,vtk_ihigh,vtk_jlow,vtk_jhigh
        integer :: nbr_left,nbr_right,nbr_bottom,nbr_top
        integer :: local_bounds(4)
        integer, allocatable :: all_bounds(:,:)
        integer :: global_ilow,global_ihigh,global_jlow,global_jhigh
        integer :: part
        character(len=40) :: timestamp_str,dirname,rank_str,q_name
        character(len=128) :: source
        character(len=200) :: filename

        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call info%get_local_block_bounds(ilow,ihigh,jlow,jhigh)
        call info%get_neighbouring_ranks(left=nbr_left,right=nbr_right,bottom=nbr_bottom,top=nbr_top)

        vtk_ilow = ilow
        vtk_ihigh = ihigh
        vtk_jlow = jlow
        vtk_jhigh = jhigh
        if(nbr_left/=MPI_PROC_NULL) vtk_ilow = ilow - 1
        if(nbr_bottom/=MPI_PROC_NULL) vtk_jlow = jlow - 1
        if((vtk_ilow<ilow .or. vtk_jlow<jlow) .and. num_ghost<1) then
            write(*,*) "Parallel VTK output requires at least one ghost layer."
            call MPI_Abort(MPI_COMM_WORLD,1,ierr)
        end if

        local_bounds = [vtk_ilow,vtk_ihigh,vtk_jlow,vtk_jhigh]
        allocate(all_bounds(4,size))

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Gather(local_bounds,4,MPI_INTEGER,all_bounds,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if(rank==0)then
            global_ilow = minval(all_bounds(1,:))
            global_ihigh = maxval(all_bounds(2,:))
            global_jlow = minval(all_bounds(3,:))
            global_jhigh = maxval(all_bounds(4,:))

            write(timestamp_str,'(F20.10)') time
            dirname = 'output/' // trim(adjustl(timestamp_str))
            filename = trim(dirname) // '/solution.pvts'

            error = pvts%initialize(filename=trim(filename),mesh_topology='PStructuredGrid',mesh_kind='Float64', &
                nx1=global_ilow,nx2=global_ihigh,ny1=global_jlow,ny2=global_jhigh,nz1=1,nz2=1)
            if(error/=0) call abort_pvts_write('Could not initialize parallel VTK file.',error)

            error = pvts%xml_writer%write_dataarray(location='node',action='open')
            if(error/=0) call abort_pvts_write('Could not open parallel VTK point data.',error)
            do part=1,num_equations
                write(q_name,'(A,I0)') 'q',part
                error = pvts%xml_writer%write_parallel_dataarray(data_name=trim(q_name),data_type='Float64', &
                    number_of_components=1_I4P)
                if(error/=0) call abort_pvts_write('Could not declare parallel VTK data array.',error)
            end do
            error = pvts%xml_writer%write_dataarray(location='node',action='close')
            if(error/=0) call abort_pvts_write('Could not close parallel VTK point data.',error)

            do part=1,size
                write(rank_str,'(I0)') part-1
                source = 'rank_' // trim(rank_str) // '/solution_rank_' // trim(rank_str) // '.vts'
                error = pvts%xml_writer%write_parallel_geo(source=trim(source), &
                    nx1=all_bounds(1,part),nx2=all_bounds(2,part), &
                    ny1=all_bounds(3,part),ny2=all_bounds(4,part),nz1=1,nz2=1)
                if(error/=0) call abort_pvts_write('Could not add a VTK rank piece.',error)
            end do

            error = pvts%finalize()
            if(error/=0) call abort_pvts_write('Could not finalize parallel VTK file.',error)

            call write_pvd_collection(time,error)
            if(error/=0) call abort_pvts_write('Could not update VTK time-series collection.',error)
        end if

        deallocate(all_bounds)
#endif
#ifdef USE_VTK
    contains
        subroutine abort_pvts_write(message,error_code)
            character(len=*), intent(in) :: message
            integer, intent(in) :: error_code
            integer :: abort_error

            write(*,*) trim(message), ' Error code: ', error_code
            call MPI_Abort(MPI_COMM_WORLD,error_code,abort_error)
        end subroutine abort_pvts_write
#endif
    end subroutine write_pvts_rank

#ifdef USE_VTK
    subroutine write_pvd_collection(time,error)
        double precision, intent(in) :: time
        integer, intent(out) :: error

        double precision, allocatable :: updated_times(:)
        integer :: count,unit,close_error,index
        character(len=256) :: message
        character(len=40) :: time_str,directory

        error = 0
        count = 0
        if(allocated(pvd_times)) count = size(pvd_times)
        if(count>0) then
            if(pvd_times(count)==time) return
        end if

        allocate(updated_times(count+1),stat=error,errmsg=message)
        if(error/=0) then
            write(*,*) 'Could not allocate VTK time-series entries: ',trim(message)
            return
        end if
        if(count>0) updated_times(1:count) = pvd_times
        updated_times(count+1) = time
        call move_alloc(updated_times,pvd_times)

        open(newunit=unit,file='output/solution.pvd',status='replace',action='write', &
            iostat=error,iomsg=message)
        if(error/=0) then
            write(*,*) 'Could not open output/solution.pvd: ',trim(message)
            return
        end if

        write(unit,'(A)',iostat=error,iomsg=message) '<?xml version="1.0"?>'
        if(error==0) write(unit,'(A)',iostat=error,iomsg=message) &
            '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">'
        if(error==0) write(unit,'(A)',iostat=error,iomsg=message) '  <Collection>'
        do index=1,count+1
            if(error/=0) exit
            write(time_str,'(F20.10)') pvd_times(index)
            directory = adjustl(time_str)
            write(unit,'(A)',iostat=error,iomsg=message) &
                '    <DataSet timestep="'//trim(directory)//'" group="" part="0" file="'// &
                trim(directory)//'/solution.pvts"/>'
        end do
        if(error==0) write(unit,'(A)',iostat=error,iomsg=message) '  </Collection>'
        if(error==0) write(unit,'(A)',iostat=error,iomsg=message) '</VTKFile>'
        close(unit,iostat=close_error,iomsg=message)
        if(error==0 .and. close_error/=0) error=close_error
        if(error/=0) write(*,*) 'Could not write output/solution.pvd: ',trim(message)
    end subroutine write_pvd_collection
#endif

    subroutine writeData(info,time,solution)
        !! Writes the solution data (q, aux) to disk.
        use mpi_f08
        use mpidcl,only:decomp_info
        use MOD_domain
        use MOD_solver_parameters

        type(decomp_info), intent(in) :: info
        double precision, intent(in) :: time
        type(solution_fields), intent(in) :: solution

        integer :: i,j,m
        character(len=40) :: filename_q
        character(len=40) :: filename_a

        character(len=40) :: timestamp_str
        character(len=40) :: dirname,subdirname,rank_str

        integer :: rank,size,ierr

        integer :: ilow,ihigh, jlow,jhigh

        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        call info%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        write(rank_str, '(I0)') rank


        ! Format to string with up to 10 digits after decimal
        write(timestamp_str, '(F20.10)') time

        ! Create parent directory name
        dirname = 'output/'//trim(adjustl(timestamp_str))
        call create_directory(dirname, rank)

        ! Create per-rank subdirectory
        subdirname = adjustl(trim(dirname) // '/rank_' // trim(rank_str))
        call create_directory(subdirname, rank)

        !create the filenames
        filename_q='/qout.dat'
        filename_a='/aout.dat'

1001    format(50e26.16)

        !give each file write operation to a separate thread (if possible)
        !$OMP PARALLEL SECTIONS
        !$OMP SECTION
        !write solution
        open(unit=1,file=adjustl(trim(subdirname))//adjustl(trim(filename_q)),status='unknown',form='formatted') !rewrite the file if exists

        do j=jlow-num_ghost,jhigh+num_ghost
            do i=ilow-num_ghost,ihigh+num_ghost
                do m=1,num_equations
                    !call padding_q(m,i,j)
                end do
                write(1,1001) (solution%vars(m,i,j), m=1,num_equations)
            end do
        end do

        close(unit=1)

        !$OMP SECTION
        !write aux
        open(unit=2,file=adjustl(trim(subdirname))//adjustl(trim(filename_a)),status='unknown',form='formatted') !rewrite the file if exists

        do j=jlow-num_ghost,jhigh+num_ghost
            do i=ilow-num_ghost,ihigh+num_ghost
                do m=1,num_aux
                    !call padding_aux(m,i,j)
                end do
                write(2,1001) (solution%aux(m,i,j), m=1,num_aux)
            end do
        end do
        close(unit=2)
        !$OMP END PARALLEL SECTIONS

    end subroutine

    subroutine write_mesh_to_disk(info,mesh)
        !! Writes the phsical coordinates of the mesh to disk.
        use mpi_f08
        use mpidcl,only:decomp_info
        use MOD_solver_parameters
        use MOD_domain

        type(decomp_info), intent(in) :: info
        type(mesh_fields), intent(in) :: mesh

        integer :: rank,ierr
        integer :: ilow,ihigh,jlow,jhigh
        character(len=40) :: rank_str
        integer :: i,j
        double precision :: x_local, y_local
        double precision :: xi_local, eta_local

        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        call info%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        write(rank_str, '(I0)') rank
        ! Write each rank's mesh to disk
        open(unit=1,file=adjustl('output/mesh_rank_' // trim(rank_str)//'.dat'),status='unknown',form='formatted')
1001    format(50e26.16)

        do j=jlow-num_ghost,jhigh+num_ghost
            do i=ilow-num_ghost,ihigh+num_ghost
                x_local = mesh%physical_space%x(i,j)
                y_local = mesh%physical_space%y(i,j)
                xi_local = mesh%computational_space%xi(i,j)
                eta_local = mesh%computational_space%eta(i,j)
                !as a safety measure: padd the values of the coordinates in case the have a extremely small value.
                !Otherwise this might result in problems when reading into matlab or similar programs.
                if (abs(x_local)<1d-99) x_local = 0.d0
                if (abs(y_local)<1d-99) y_local = 0.d0
                if (abs(xi_local)<1d-99) xi_local = 0.d0
                if (abs(eta_local)<1d-99) eta_local = 0.d0
                write(1,1001) mesh%physical_space%x(i,j), mesh%physical_space%y(i,j), mesh%computational_space%xi(i,j), mesh%computational_space%eta(i,j)
            end do
        end do
        close(unit=1)
    end subroutine

    subroutine create_directory(dir, rank)
        !! Helper to create a directory.
        use mpi_f08

        character(len=*), intent(in) :: dir
        integer, intent(in) :: rank

        integer :: stat,ierr

        call execute_command_line('mkdir -p ' // trim(dir), exitstat=stat)
        if (stat /= 0) then
            print *, 'Rank ',rank,' : Failed to create directory:', trim(dir)
        end if

        ! Synchronize to ensure directories are created before use
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine create_directory

end module MOD_IO
