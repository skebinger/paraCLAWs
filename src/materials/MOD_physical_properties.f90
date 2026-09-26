module MOD_materials
    implicit none

    type physical_property
        double precision, allocatable :: value (:,:,:)
    contains
        procedure :: set_val
        procedure :: get_val
    end type physical_property

    type material
        character(len=128) :: name
        type(physical_property), allocatable :: properties(:)
    end type material

    interface allocate
        module procedure allocate_property
        module procedure allocate_material
    end interface allocate

    interface
        module subroutine allocate_property(this, decomposition)
            use mpidcl,only: decomp_info
            type(physical_property), intent(inout) :: this
            type(decomp_info), intent(in) :: decomposition
        end subroutine allocate_property

        module subroutine allocate_material(this, decomposition, num_materials)
            use mpidcl,only: decomp_info
            type(material), intent(inout) :: this
            type(decomp_info), intent(in) :: decomposition
            integer, intent(in) :: num_materials
        end subroutine allocate_material
    end interface

    interface
        module subroutine set_val(this, val)
            class(physical_property), intent(inout) :: this
            double precision, intent(in) :: val
        end subroutine set_val

        module function get_val(this) result(val)
            class(physical_property), intent(in) :: this
            double precision :: val
        end function get_val
    end interface


    type(material), allocatable :: materials(:)

contains



end module MOD_materials
