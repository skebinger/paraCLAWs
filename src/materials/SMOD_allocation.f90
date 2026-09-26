submodule(MOD_materials) SMOD_allocation
    implicit none

contains

    module subroutine allocate_property(property, decomposition)
        use mpidcl,only: decomp_info
        type(physical_property), intent(inout) :: property
        type(decomp_info), intent(in) :: decomposition

        integer :: ilow,ihigh,jlow,jhigh,klow,khigh

        call decomposition%get_local_block_bounds(ilow,ihigh,jlow,jhigh,klow,khigh)

        allocate(property%value(ilow:ihigh, jlow:jhigh, klow:khigh))
    end subroutine allocate_property

    module subroutine allocate_material(material, decomposition, num_materials)
        use mpidcl,only: decomp_info
        type(material), intent(inout) :: material
        type(decomp_info), intent(in) :: decomposition
        integer, intent(in) :: num_materials

        ! Allocate the properties array for the material
        ! Assuming some number of properties, e.g., 3 for demonstration
        allocate(material%properties(num_materials))

        ! Allocate each property
        ! Note: You might need to pass decomposition to each property's allocate
        ! But since allocate_property takes decomposition, we can call it
        ! However, in this context, material might need to handle its own allocation logic
    end subroutine allocate_material

end submodule
