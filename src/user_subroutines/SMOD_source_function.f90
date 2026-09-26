! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_user_subroutines) SMOD_source_function
    implicit none

contains

    module function source_function(i,j,mesh,solution,decomposition)result(phi)
        use MOD_user_data
        use MOD_solver_parameters,only:num_equations
        use MOD_index
        use MOD_domain,only:mesh_fields,solution_fields
        use mpidcl,only:decomp_info
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(in) :: solution
        type(decomp_info), intent(in) :: decomposition
        integer :: i,j
        double precision :: phi(num_equations)

        phi = 0.0d0
    end function


end submodule
