! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_flux_calculation) SMOD_approximate_Riemann_solvers
    implicit none

contains

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

        !================================================================================
        !HLL solver variables (2-wave solver); resolving contact discontinuity of equations
        double precision :: q_intermediate(num_equations)
        double precision :: flux_function_left(num_equations)
        double precision :: flux_function_right(num_equations)
        double precision :: jakobi_left(num_equations,num_equations) !! jakobi of left state (i-1)
        double precision :: jakobi_right(num_equations,num_equations) !! jakobi of right state (i)
        double precision :: lambda_left(num_equations) !! eigenvalues of left state (i-1)
        double precision :: lambda_right(num_equations) !! eigenvalues of right state (i)
        double precision :: s1 !! wave speed 1; named for consistency with LeVeque FVMHP book
        double precision :: s2 !! wave speed 2; named for consistency with LeVeque FVMHP book
        double precision :: R_left(num_equations,num_equations) !! matrix of eigenvectors R=[r_1,r_2,...,r_n] of left state (i-1)
        double precision :: R_right(num_equations,num_equations) !! matrix of eigenvectors R=[r_1,r_2,...,r_n] of right state (i)

        if(ixy==1)then
            !get the left state jakobi (i-1)
            jakobi_left=getJakobi(ixy,i-1,j,mesh,solution,work,decomposition)

            !get the right state jakobi (i)
            jakobi_right=getJakobi(ixy,i,j,mesh,solution,work,decomposition)
        else
            !get the left state jakobi (j-1)
            jakobi_left=getJakobi(ixy,i,j-1,mesh,solution,work,decomposition)

            !get the right state jakobi (j)
            jakobi_right=getJakobi(ixy,i,j,mesh,solution,work,decomposition)
        end if

        !now calculate the eigenvalues, eigenvectors for these two states
        !this can be made faster, if J_left is initialized at the beginning and at the end J_right is written into that
        !-> allows for skipping one call to this subroutine
        call calculateEigenSolution(num_equations,num_equations,jakobi_left,R_left,lambda_left,.True.)
        call calculateEigenSolution(num_equations,num_equations,jakobi_right,R_right,lambda_right,.True.)

        !find the biggest eigenvalue
        s2=maxval([abs(lambda_left),abs(lambda_right)])
        work%lambda(2,i,j) = s2

        !find the smalles eigenvalue
        s1=-s2!minval([lambda_left,lambda_right])
        work%lambda(1,i,j) = s1

        !initialize flux function variables before inputting them into the matrix vector product (SAFETY REASON; might cause uninit. values otherwise)
        flux_function_left=0.d0
        flux_function_right=0.d0

        !check for both eigenvalues = 0 (allows to skip this section and set both waves to 0)
        if(abs(work%lambda(1,i,j))>ud_precision.AND.abs(work%lambda(2,i,j))>ud_precision)then
            !calculate the left and right flux functions (f_i = jakobian_i * q_i, since jakobian = df_dq)
            if(ixy==1)then
                call matrix_vector_product(num_equations,num_equations,jakobi_left,solution%vars(:,i-1,j),flux_function_left)
                call matrix_vector_product(num_equations,num_equations,jakobi_right,solution%vars(:,i,j),flux_function_right)
            else
                call matrix_vector_product(num_equations,num_equations,jakobi_left,solution%vars(:,i,j-1),flux_function_left)
                call matrix_vector_product(num_equations,num_equations,jakobi_right,solution%vars(:,i,j),flux_function_right)
            end if

            !calculate the intermediate state @ i-1/2
            if(ixy==1)then
                q_intermediate = (flux_function_right - flux_function_left - s2*solution%vars(:,i,j)  + s1*solution%vars(:,i-1,j)) / (ud_precision + s1-s2)
                work%waves(:,1,i,j) = (q_intermediate - solution%vars(:,i-1,j))
            else
                q_intermediate = (flux_function_right - flux_function_left - s2*solution%vars(:,i,j)  + s1*solution%vars(:,i,j-1)) / (ud_precision + s1-s2)
                work%waves(:,1,i,j) = (q_intermediate - solution%vars(:,i,j-1))
            end if
            work%waves(:,2,i,j) = (solution%vars(:,i,j) - q_intermediate)
        else
            work%lambda(:,i,j) = 0.d0
            work%waves(:,:,i,j) = 0.d0
        end if
    end subroutine

    !module subroutine ROE_like()
    !end subroutine


end submodule
