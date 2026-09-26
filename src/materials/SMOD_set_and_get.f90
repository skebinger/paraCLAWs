submodule(MOD_materials) SMOD_set_and_get
    implicit none

contains

    module subroutine set_val(this, val)
        class(physical_property), intent(inout) :: this
        double precision, intent(in) :: val

        this%value = val
    end subroutine set_val

    module function get_val(this) result(val)
        class(physical_property), intent(in) :: this
        double precision :: val
        
        val = this%value
    end function get_val

end submodule
