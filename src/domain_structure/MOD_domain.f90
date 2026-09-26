! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_domain
    !! Module containing all field variables and domain information.
    implicit none

    !===========================================================
    ! Types
    !===========================================================
    type mesh_base
        !! Base mesh information, independent of the actual solution.
        integer :: m_xi !! number of cells in xi direction (in computational space)
        integer :: m_eta !! number of cells in eta direction (in computational space)
        double precision :: xi_dimensions(2) !! xi dimensions of the domain (in computational space)
        double precision :: eta_dimensions(2) !! eta dimensions of the domain (in computational space)

    contains
        procedure :: setup_mesh_base
    end type mesh_base

    !===========================================================
    type mesh_face
        !! so far an empty unused type; might be picked up if implicit solver is added in the future
    end type

    !===========================================================
    type, extends(mesh_face) :: mesh_fields_computational_space
        !! Mesh information in computational space, including the actual coordinates and geometric factors.
        double precision :: d_xi !! equidistant spacing in xi direction (in computational space)
        double precision :: d_eta !! equidistant spacing in eta direction (in computational space)

        double precision, allocatable :: xi(:,:) !! xi coordinate (in computational space)
        double precision, allocatable :: eta(:,:) !! eta coordinate (in computational space)
        double precision, allocatable :: xi_corner(:,:,:) !! xi corner coordinate (in computational space)
        double precision, allocatable :: eta_corner(:,:,:) !! eta corner coordinate (in computational space)

    contains
        procedure :: allocate_computational_space
        procedure :: deallocate_computational_space
        procedure :: initialize_computational_coordinates
    end type

    !===========================================================
    type, extends(mesh_face) :: mesh_fields_physical_space
        !! Mesh information in physical space, including the actual coordinates and geometric factors.
        double precision, allocatable :: x(:,:) !! x coordinate (in physical space)
        double precision, allocatable :: y(:,:) !! y coordinate (in physical space)
        double precision, allocatable :: x_corner(:,:,:) !! x corner coordinate (in physical space)
        double precision, allocatable :: y_corner(:,:,:) !! y corner coordinate (in physical space)

    contains
        procedure :: allocate_physical_space
        procedure :: deallocate_physical_space
    end type

    !===========================================================
    type mesh_fields_quadrilateral_mapping
        !! Geometric factors of the mesh, needed for the solver.
        double precision, allocatable :: capacity(:,:) !! capacity array
        double precision, allocatable :: normal_vector_xi(:,:,:,:) !! normal vector of a xi-edge
        double precision, allocatable :: normal_vector_eta(:,:,:,:) !! normal vector of an eta-edge
        double precision, allocatable :: length_ratio_xi(:,:,:) !! length ratio of an edge along xi (dy/deta)
        double precision, allocatable :: length_ratio_eta(:,:,:) !! length ratio of an edge along eta (dx/dxi)
    contains
        procedure :: allocate_quadrilateral_mapping
        procedure :: deallocate_quadrilateral_mapping
    end type

    !===========================================================
    type, extends(mesh_base) :: mesh_fields
        !! Type containing all mesh information, including the actual coordinates and geometric factors.
        type(mesh_fields_computational_space) :: computational_space
        type(mesh_fields_physical_space) :: physical_space
        type(mesh_fields_quadrilateral_mapping) :: quadrilateral_mapping
    contains
        procedure :: finalize_mesh
    end type mesh_fields

    !===========================================================
    type solution_fields
        !! Type containing all solution information, including the actual solution and auxiliary fields.
        double precision, allocatable :: vars(:,:,:) !! The solution variables
        double precision, allocatable :: aux(:,:,:) !! Auxiliary variables
    contains
        procedure :: allocate_solution
        procedure :: deallocate_solution
    end type solution_fields

    !===========================================================
    ! Type bound procedure interfaces
    !===========================================================
    ! Allocation interfaces
    interface
        module subroutine allocate_computational_space(this,decomposition)
            use mpidcl,only:decomp_info
            use MOD_solver_parameters, only: num_ghost
            class(mesh_fields_computational_space), intent(inout) :: this
            type(decomp_info), intent(in) :: decomposition
        end subroutine

        module subroutine allocate_physical_space(this,decomposition)
            use mpidcl,only:decomp_info
            use MOD_solver_parameters, only: num_ghost
            class(mesh_fields_physical_space), intent(inout) :: this
            type(decomp_info), intent(in) :: decomposition
        end subroutine

        module subroutine allocate_quadrilateral_mapping(this,decomposition)
            use mpidcl,only:decomp_info
            use MOD_solver_parameters, only: num_ghost
            class(mesh_fields_quadrilateral_mapping), intent(inout) :: this
            type(decomp_info), intent(in) :: decomposition
        end subroutine

        module subroutine allocate_solution(this,decomposition)
            use mpidcl,only:decomp_info
            use MOD_solver_parameters,only:num_equations,num_ghost,num_aux
            class(solution_fields), intent(inout) :: this
            type(decomp_info), intent(in) :: decomposition
        end subroutine
    end interface

    ! Deallocation interfaces
    interface
        module subroutine deallocate_computational_space(this)
            class(mesh_fields_computational_space), intent(inout) :: this
        end subroutine

        module subroutine deallocate_physical_space(this)
            class(mesh_fields_physical_space), intent(inout) :: this
        end subroutine

        module subroutine deallocate_quadrilateral_mapping(this)
            class(mesh_fields_quadrilateral_mapping), intent(inout) :: this
        end subroutine

        module subroutine deallocate_solution(this)
            class(solution_fields), intent(inout) :: this
        end subroutine
    end interface


    ! Interface for mesh initialization and manipulation
    interface
        module subroutine initialize_computational_coordinates(c_mesh,decomposition,m_xi,m_eta,xi_dimensions,eta_dimensions)
            !! Initialize the computational mesh coordinates at the cell centre and cell corners.
            use MOD_solver_parameters,only:num_ghost
            use mpidcl,only:decomp_info
            class(mesh_fields_computational_space), intent(inout) :: c_mesh !! the mesh in computational space
            type(decomp_info), intent(in) :: decomposition
            integer, intent(in) :: m_xi,m_eta
            double precision, intent(in) :: xi_dimensions(2), eta_dimensions(2)
        end subroutine

        module subroutine finalize_mesh(mesh,decomposition)
            !! Finalize the mesh coordinate fields and the quadrilateral mapping information.
            !!
            !! First fill the computational coordinates and second map
            !! computational to physical coordinates (TBD).
            use mpidcl,only:decomp_info
            use MOD_solver_parameters,only:num_ghost

            class(mesh_fields), intent(inout) :: mesh
            type(decomp_info), intent(in) :: decomposition
        end subroutine

        module subroutine rescale_coordinates(mesh,decomposition)
            use mpidcl,only:decomp_info
            use MOD_user_data
            class(mesh_fields), intent(inout) :: mesh
            type(decomp_info), intent(in) :: decomposition
        end subroutine

        module subroutine map_to_physical_coordinates(xp,yp,xi,eta)
            double precision, intent(out) :: xp,yp
            double precision, intent(in) :: xi,eta
        end subroutine

    end interface
    interface
        module subroutine update_rank_boundaries(decomposition,solution)
            use mpi_f08
            use mpidcl,only:exchange_halos,get_cartesian_comm,get_neighbouring_ranks,decomp_info
            use MOD_solver_parameters,only:num_equations,num_ghost,num_aux

            type(decomp_info), intent(in) :: decomposition
            type(solution_fields), intent(inout) :: solution
        end subroutine
    end interface

contains

    !================================================================================================

    subroutine setup_mesh_base(this,m_xi,m_eta,xi_dimensions,eta_dimensions)
        class(mesh_base), intent(inout) :: this
        integer, intent(in) :: m_xi,m_eta
        double precision, intent(in) :: xi_dimensions(2),eta_dimensions(2)

        this%m_xi=m_xi
        this%m_eta=m_eta

        this%xi_dimensions=xi_dimensions
        this%eta_dimensions=eta_dimensions
    end subroutine

    subroutine write_base_info(mesh)
        use mpi_f08
        use MOD_solver_parameters,only:num_ghost,num_aux,num_equations
        type(mesh_fields), intent(in) :: mesh

        integer :: rank,size,ierr
        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

        ! write base info to disk, needed for post processing the data
        if(rank==0)then
            open(unit=1,file="./output/domain.dat",status='unknown',form='formatted')
            write(1,*) size
            write(1,*) mesh%m_xi
            write(1,*) mesh%m_eta
            write(1,*) num_equations
            write(1,*) num_aux
            write(1,*) num_ghost
            write(1,*) mesh%xi_dimensions(1)
            write(1,*) mesh%xi_dimensions(2)
            write(1,*) mesh%eta_dimensions(1)
            write(1,*) mesh%eta_dimensions(2)
            close(unit=1)
        end if

    end subroutine

end module MOD_domain
