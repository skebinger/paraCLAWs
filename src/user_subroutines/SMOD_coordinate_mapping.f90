! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_domain) SMOD_coordinate_mapping
    implicit none

contains

    module subroutine rescale_coordinates(mesh,decomposition)
        !! Rescale the computational coordinates.
        !! This is thought to be a method, for setting the correct final computational mesh prior to calculating the quadrilateral
        !! mapping information if specifying the correct mesh bounds (lower and upper bounds) in the control.in file is not easily done.
        !! For example: You know it is a multiple of some quantity you calculate inside the user setup subroutine, but don't know its
        !! value at compile time.
        use MOD_solver_parameters,only:num_ghost
        use mpidcl,only:decomp_info
        use MOD_user_data
        class(mesh_fields), intent(inout) :: mesh
        type(decomp_info), intent(in) :: decomposition

        integer :: ilow, ihigh, jlow, jhigh
        integer :: i,j

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh)

    end subroutine

    module subroutine map_to_physical_coordinates(xp,yp,xi,eta)
        !! Subroutine to map from computational coordinates to physical coordinates.
        double precision, intent(out) :: xp,yp
        double precision, intent(in) :: xi,eta

        xp = xi
        yp = eta

    end subroutine

end submodule
