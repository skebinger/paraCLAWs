! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_flux_calculation
    !! Module contains all neccessary routines to calculate the approximate Riemann solution
    !! for a given cell by decomposing the jump at the interface i-1/2 or j-1/2 into waves
    !! and according wave speeds.
    implicit none

    INTERFACE
        module function getJakobi(ixy,i,j,mesh,solution,work,decomposition)result(jakobi)
            !! Returns the jakobian matrix for the state vars(i,j); w/ aux array
            use mpi_f08
            use MOD_solver_parameters,only:num_equations,psi=>grav
            use MOD_user_subroutines
            use MOD_user_data
            use MOD_index
            use MOD_domain,only:solution_fields,mesh_fields
            use MOD_mathematical_functions, only: integrate
            use MOD_work_arrays,only:work_arrays
            use mpidcl,only:decomp_info
            integer, intent(in) :: ixy !! direction of sweep
            integer, intent(in) :: i,j
            type(mesh_fields), intent(in) :: mesh
            type(solution_fields), intent(in) :: solution
            type(work_arrays), intent(in) :: work
            type(decomp_info), intent(in) :: decomposition
            double precision :: jakobi(num_equations,num_equations)
        end function
    END INTERFACE

    interface
        module subroutine HLL(ixy,i,j,mesh,solution,work,decomposition)
            !! Executes the calculation of waves and wave speeds using the HLL approximate Riemann solver.
            use MOD_solver_parameters,only:num_equations
            use MOD_constants,only:ud_precision
            use MOD_mathematical_functions,only:inverse_matrix,matrix_vector_product,calculateEigenSolution
            use MOD_domain,only:mesh_fields,solution_fields
            use MOD_work_arrays,only:work_arrays
            use mpidcl,only:decomp_info

            !================================================================================
            integer, intent(in) :: ixy !! direction of current sweep
            integer, intent(in) :: i,j
            type(mesh_fields), intent(in) :: mesh
            type(solution_fields), intent(inout) :: solution
            type(work_arrays), intent(inout) :: work
            type(decomp_info), intent(in) :: decomposition
            !================================================================================
        end subroutine
    end interface

contains

    subroutine wave_decomposition(ixy,i,j,mesh,solution,work,decomposition)
        !! Subroutine executing the decomposition. Acts mostly as a wrapper for the chosen
        !! Riemann solver.
        use MOD_domain,only:mesh_fields,solution_fields
        use MOD_work_arrays,only:work_arrays
        use mpidcl,only:decomp_info
        !================================================================================
        integer, intent(in) :: ixy !! direction of current sweep
        integer, intent(in) :: i,j
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution
        type(work_arrays), intent(inout) :: work
        type(decomp_info), intent(in) :: decomposition
        !================================================================================


        call HLL(ixy,i,j,mesh,solution,work,decomposition)

    end subroutine


    subroutine split_flux(decomposition,lambda,waves,Aminus_dQ,Aplus_dQ)
        !! Split the waves into left and right going waves.
        !! The contributions are split into the components Aminus_dQ and Aplus_dQ,
        !! which are used to update the solution in the dimensional splitting step. <br>
        !! $$A^- \Delta Q = \sum_{mw} \min(\lambda_{mw},0) W_{mw}$$ <br>
        !! $$A^+ \Delta Q = \sum_{mw} \max(\lambda_{mw},0) W_{mw}$$
        use mpidcl,only:decomp_info
        use MOD_solver_parameters,only:num_equations,num_waves

        type(decomp_info), intent(in) :: decomposition
        double precision, allocatable, intent(in) :: lambda(:,:,:) !! vector of wave speeds for each wave and cell
        double precision, allocatable, intent(in) :: waves(:,:,:,:) !! waves split into wave families for each cell
        double precision, allocatable, intent(inout) :: Aminus_dQ(:,:,:) !! contribution of left-going waves to the update
        double precision, allocatable, intent(inout) :: Aplus_dQ(:,:,:) !! contribution of right-going waves to the update

        integer :: i,j
        integer :: me,mw
        integer :: ilow, ihigh, jlow, jhigh

        double precision :: s,w,mask

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(mask,s,w,i,mw,me)
        do j=jlow-1,jhigh+1
            do i=ilow-1,ihigh+1
                do mw=1,num_waves
                    s=lambda(mw,i,j)
                    do me=1,num_equations
                        w=waves(me,mw,i,j)
                        mask = merge(1.d0, 0.d0, s<0.d0)

                        Aminus_dQ(me,i,j) = Aminus_dQ(me,i,j) + mask*s*w
                        Aplus_dQ(me,i,j) = Aplus_dQ(me,i,j) +(1.d0-mask)*s*w
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine

end module
