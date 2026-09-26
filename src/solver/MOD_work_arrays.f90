! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_work_arrays
    use MOD_look_up_table,only:lu_tab_1D
    implicit none

    type work_arrays
        !! A collection of work arrays that are constantly used
        double precision, allocatable :: lambda(:,:,:)
        double precision, allocatable :: waves(:,:,:,:)
        double precision, allocatable :: Aminus_dQ(:,:,:)
        double precision, allocatable :: Aplus_dQ(:,:,:)

        double precision, allocatable :: wnorm2l(:,:,:)
        double precision, allocatable :: wnorm2r(:,:,:)
        double precision, allocatable :: f(:,:,:)

        type(lu_tab_1D) :: beta_tab_func
        type(lu_tab_1D) :: beta_tab_deriv
        type(lu_tab_1D) :: F_tab
        type(lu_tab_1D) :: G_tab
        type(lu_tab_1D) :: dF_tab
        type(lu_tab_1D) :: dG_tab
    contains
        procedure :: allocate_work_arrays
        procedure :: deallocate_work_arrays
        procedure :: flash_work_to_zero
        procedure :: initialize_work
    end type

    interface
        module subroutine initialize_work(work,decomposition)
            use mpidcl,only:decomp_info
            class(work_arrays) :: work
            type(decomp_info), intent(in) :: decomposition
        end subroutine
    end interface

contains

    subroutine allocate_work_arrays(work,info)
        use MOD_solver_parameters,only:num_waves,num_equations,num_ghost
        use mpidcl,only:get_local_block_bounds,decomp_info

        class(work_arrays), intent(out) :: work
        type(decomp_info), intent(in) :: info

        integer :: xi_start, xi_end, eta_start, eta_end

        ! get MPI block bounds
        call info%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        ! first order integration loop variables
        allocate(work%lambda(num_waves,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))
        allocate(work%waves(num_equations,num_waves,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))
        allocate(work%Aminus_dQ(num_equations,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))
        allocate(work%Aplus_dQ(num_equations,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))

        ! higher order correction
        allocate(work%f(num_equations,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))
        allocate(work%wnorm2l(num_waves,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))
        allocate(work%wnorm2r(num_waves,xi_start-num_ghost:xi_end+num_ghost,eta_start-num_ghost:eta_end+num_ghost))
    end subroutine

    subroutine deallocate_work_arrays(work)
        class(work_arrays), intent(inout) :: work

        deallocate(work%lambda,work%waves,work%Aminus_dQ,work%Aplus_dQ,work%f,work%wnorm2l,work%wnorm2r)
    end subroutine

    subroutine flash_work_to_zero(work,decomposition)
        use MOD_solver_parameters,only:num_ghost,num_equations,num_waves
        use mpidcl,only:decomp_info
        type(decomp_info), intent(in) :: decomposition
        class(work_arrays), intent(inout) :: work

        integer :: i,j,me,mw
        integer :: xi_start, xi_end, eta_start, eta_end

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        !$OMP PARALLEL DO PRIVATE(i,mw,me) COLLAPSE(2)
        do j=eta_start-num_ghost,eta_end+num_ghost
            do i=xi_start-num_ghost,xi_end+num_ghost
                do mw=1,num_waves

                    work%wnorm2l(mw,i,j)=0.d0
                    work%wnorm2r(mw,i,j)=0.d0
                    work%lambda(mw,i,j)=0.d0

                    do me=1,num_equations
                        work%Aminus_dQ(me,i,j)=0.d0
                        work%Aplus_dQ(me,i,j)=0.d0
                        work%waves(me,mw,i,j)=0.d0
                        work%f(me,i,j)=0.d0
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine
end module
