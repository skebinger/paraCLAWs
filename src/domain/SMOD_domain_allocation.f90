! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_domain) SMOD_domain_allocation
    implicit none

contains

    module subroutine allocate_computational_space(this,decomposition)
        use mpidcl,only:decomp_info
        use MOD_solver_parameters, only: num_ghost
        class(mesh_fields_computational_space), intent(inout) :: this
        type(decomp_info), intent(in) :: decomposition

        integer :: ilow, ihigh, jlow, jhigh

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        allocate(this%xi (ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%eta(ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%xi_corner (5,ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%eta_corner(5,ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))

        ! Initialize
        this%xi(:,:) = 0
        this%eta(:,:) = 0
        this%xi_corner(:,:,:) = 0
        this%eta_corner(:,:,:) = 0
    end subroutine

    module subroutine allocate_physical_space(this,decomposition)
        use mpidcl,only:decomp_info
        use MOD_solver_parameters, only: num_ghost
        class(mesh_fields_physical_space), intent(inout) :: this
        type(decomp_info), intent(in) :: decomposition

        integer :: ilow, ihigh, jlow, jhigh

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        allocate(this%x(ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%y(ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%x_corner(5,ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%y_corner(5,ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))

        ! Initialize
        this%x(:,:) = 0
        this%y(:,:) = 0
        this%x_corner(:,:,:) = 0
        this%y_corner(:,:,:) = 0
    end subroutine

    module subroutine allocate_quadrilateral_mapping(this,decomposition)
        use mpidcl,only:decomp_info
        use MOD_solver_parameters, only: num_ghost
        class(mesh_fields_quadrilateral_mapping), intent(inout) :: this
        type(decomp_info), intent(in) :: decomposition

        integer :: ilow, ihigh, jlow, jhigh

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        allocate(this%capacity          (      ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%normal_vector_xi  (2, 2, ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%normal_vector_eta (2, 2, ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%length_ratio_xi   (2, ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%length_ratio_eta  (2, ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))

        ! Initialize
        this%capacity(:,:) = 0
        this%normal_vector_xi(:,:,:,:) = 0
        this%normal_vector_eta(:,:,:,:) = 0
        this%length_ratio_xi(:,:,:) = 0
        this%length_ratio_eta(:,:,:) = 0

    end subroutine

    module subroutine allocate_solution(this,decomposition)
        use mpidcl,only:decomp_info
        use MOD_solver_parameters,only:num_equations,num_ghost,num_aux
        class(solution_fields), intent(inout) :: this
        type(decomp_info), intent(in) :: decomposition

        integer :: ilow, ihigh, jlow, jhigh

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        allocate(this%vars  (num_equations  ,ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))
        allocate(this%aux   (num_aux        ,ilow-num_ghost:ihigh+num_ghost,jlow-num_ghost:jhigh+num_ghost))

        ! Initialize
        this%vars   = 0
        this%aux    = 0
    end subroutine

    !TBD: Deallocation subroutines

    module subroutine deallocate_computational_space(this)
        class(mesh_fields_computational_space), intent(inout) :: this
        deallocate(this%xi,this%eta,this%xi_corner,this%eta_corner)
    end subroutine

    module subroutine deallocate_physical_space(this)
        class(mesh_fields_physical_space), intent(inout) :: this
        deallocate(this%x,this%y,this%x_corner,this%y_corner)
    end subroutine

    module subroutine deallocate_quadrilateral_mapping(this)
        class(mesh_fields_quadrilateral_mapping), intent(inout) :: this
        deallocate(this%capacity,this%normal_vector_xi,this%normal_vector_eta,this%length_ratio_xi,this%length_ratio_eta)
    end subroutine

    module subroutine deallocate_solution(this)
        class(solution_fields), intent(inout) :: this
        deallocate(this%vars,this%aux)
    end subroutine


end submodule
