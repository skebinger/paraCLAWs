! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !


submodule(MOD_user_subroutines) SMOD_user_post
    implicit none
contains


    module subroutine user_post(i,j,mesh,solution,decomposition)
        !! Executes post processing at the cell (i,j)
        use MOD_index
        use MOD_domain,only:mesh_fields,solution_fields
        use mpidcl,only:decomp_info
        integer, intent(in) :: i,j
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(inout) :: solution
        type(decomp_info), intent(in) :: decomposition

        print *, "NOTHING DONE IN user_post: User must define any post processing operations to be performed at each cell (i,j)."
    end subroutine


end submodule
