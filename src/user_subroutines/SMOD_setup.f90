! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_user_subroutines) SMOD_setup

    implicit none

contains
    module subroutine setup_user_input()
        use MOD_solver_parameters,only:grav
        use MOD_constants,only:PI
        use MOD_user_data
        use MOD_mathematical_functions

        grav = 9.81

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

        integer :: ilow, ihigh, jlow, jhigh

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

        print *, "NOTHING DONE IN after_step: User must define any operations to be performed after each time step."

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

        integer :: xi_start, xi_end, eta_start, eta_end

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

    end subroutine

end submodule
