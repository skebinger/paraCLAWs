! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_user_subroutines
    implicit none

    interface
        module subroutine setup_user_input()
            use MOD_solver_parameters
            use MOD_constants
            use MOD_user_data
        end subroutine

        module subroutine after_step(decomposition,t_ctrl,mesh,solution)
            ! This subroutine can be used to perform any necessary operations after each time step.
            use mpidcl,only:decomp_info
            use MOD_time_control,only:time_control
            use MOD_domain,only:solution_fields,mesh_fields
            use MOD_index
            use MOD_user_data
            type(decomp_info), intent(in) :: decomposition
            type(time_control), intent(in) :: t_ctrl
            type(mesh_fields), intent(inout) :: mesh
            type(solution_fields), intent(inout) :: solution
        end subroutine

        module subroutine setup_auxiliary_quantities(decomposition,mesh,solution)
            use mpidcl,only:decomp_info
            use MOD_domain,only:solution_fields,mesh_fields
            use MOD_solver_parameters,only:num_ghost
            use MOD_user_data
            use MOD_index

            type(decomp_info), intent(in) :: decomposition
            type(mesh_fields), intent(in) :: mesh
            type(solution_fields), intent(inout) :: solution
        end subroutine

        module subroutine set_initial_conditions(decomposition,mesh,solution)
            use mpidcl,only:decomp_info
            use MOD_domain,only:solution_fields,mesh_fields
            use MOD_solver_parameters,only:num_ghost
            use MOD_user_data
            use MOD_index
            use mpi_f08

            type(decomp_info), intent(in) :: decomposition
            type(mesh_fields), intent(in) :: mesh
            type(solution_fields), intent(inout) :: solution
        end subroutine

        module subroutine user_post(i,j,mesh,solution,decomposition)
            !! Executes post processing at the cell (i,j)
            use MOD_index
            use MOD_domain,only:mesh_fields,solution_fields
            use mpidcl,only:decomp_info
            integer, intent(in) :: i,j
            type(mesh_fields), intent(in) :: mesh
            type(solution_fields), intent(inout) :: solution
            type(decomp_info), intent(in) :: decomposition
        end subroutine

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
        end function
    end interface

contains

end module MOD_user_subroutines
