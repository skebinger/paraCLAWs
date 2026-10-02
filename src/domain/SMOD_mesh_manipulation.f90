! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule(MOD_domain) SMOD_mesh_manipulation
    implicit none
contains

    module subroutine initialize_computational_coordinates(c_mesh,decomposition,m_xi,m_eta,xi_dimensions,eta_dimensions)
        !! Initialize the computational mesh coordinates at the cell centre and cell corners.
        use MOD_solver_parameters,only:num_ghost
        use mpidcl,only:decomp_info
        class(mesh_fields_computational_space), intent(inout) :: c_mesh !! the mesh in computational space
        type(decomp_info), intent(in) :: decomposition
        integer, intent(in) :: m_xi,m_eta
        double precision, intent(in) :: xi_dimensions(2), eta_dimensions(2)

        integer :: xi_start,xi_end,eta_start,eta_end
        integer :: i,j

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        c_mesh%d_xi = (xi_dimensions(2)-xi_dimensions(1))/m_xi
        c_mesh%d_eta = (eta_dimensions(2)-eta_dimensions(1))/m_eta

        !$OMP PARALLEL DO PRIVATE(i) COLLAPSE(2)
        do j=eta_start-num_ghost,eta_end+num_ghost
            do i=xi_start-num_ghost,xi_end+num_ghost
                ! Attention: the corners are sorted in counter-clow-wise orientation.
                ! Only then the area of the quadrilateral cell is >0!
                c_mesh%xi_corner(1,i,j) = xi_dimensions(1) + (i-1)*c_mesh%d_xi
                c_mesh%xi_corner(2,i,j) = c_mesh%xi_corner(1,i,j) + c_mesh%d_xi
                c_mesh%xi_corner(3,i,j) = c_mesh%xi_corner(2,i,j)
                c_mesh%xi_corner(4,i,j) = c_mesh%xi_corner(1,i,j)
                ! set the corner index "5" to match the point "1" -> later useful for polygon area formula:
                c_mesh%xi_corner(5,i,j) = c_mesh%xi_corner(1,i,j)

                c_mesh%eta_corner(1,i,j) = eta_dimensions(1) + (j-1)*c_mesh%d_eta
                c_mesh%eta_corner(2,i,j) = c_mesh%eta_corner(1,i,j)
                c_mesh%eta_corner(3,i,j) = c_mesh%eta_corner(2,i,j) + c_mesh%d_eta
                c_mesh%eta_corner(4,i,j) = c_mesh%eta_corner(3,i,j)
                ! set the corner index "5" to match the point "1" -> later useful for polygon area formula:
                c_mesh%eta_corner(5,i,j) = c_mesh%eta_corner(1,i,j)

                c_mesh%xi(i,j) = xi_dimensions(1) + (i-1)*c_mesh%d_xi+c_mesh%d_xi/2
                c_mesh%eta(i,j) = eta_dimensions(1) + (j-1)*c_mesh%d_eta+c_mesh%d_eta/2
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine

    module subroutine finalize_mesh(mesh,decomposition)
        !! Finalize the mesh coordinate fields and the quadrilateral mapping information.
        !!
        !! NOCH UMSCHREIBEN!!!!! EVTL MEHR MIT LATEX??
        !! This routine:
        !! 1. Maps computational coordinates to physical space via
        !!    `map_to_physical_coordinates`.
        !! 2. Computes ξ- and η-face normal vectors.
        !! 3. Computes the physical cell area using a polygon summation rule.
        !! 4. Computes the capacity:
        !!
        !!       capacity = A_physical / (Δξ · Δη)
        !!
        !! ### Corner numbering (counter-clockwise)
        !!
        !! Corners are defined as:
        !! • (1): (ξ₁,          η₁)
        !! • (2): (ξ₁ + Δξ,     η₁)
        !! • (3): (ξ₁ + Δξ,     η₁ + Δη)
        !! • (4): (ξ₁,          η₁ + Δη)
        !! • (5): = (1)  (polygon loop closure)
        !!
        !! Corresponding geometric layout (CCW ordering):
        !!
        !! ```
        !!                    η ↑
        !!
        !!            (4) ───────────── (3)
        !!             │                 │
        !!             │       ●         │    ● = cell center
        !!             │                 │
        !!            (1) ───────────── (2)
        !!                    ξ →
        !! ```
        !!
        !! ### Face vectors and outward normals
        !!
        !! The face vectors used for normal construction are:
        !! • ξ-face:  from corner 1 → 4  (vertical)
        !! • η-face:  from corner 1 → 2  (horizontal)
        !!
        !! For a face vector **v** = (dx, dy):
        !!
        !!       n = ( -dy,  dx ) / √(dx² + dy²)
        !!
        !! This yields a left-hand rotated unit vector with respect to the face.
        !!
        !! ### Polygon-based cell area
        !!
        !! The area of the cell is evaluated using the trapezoidal (shoelace) formula:
        !!
        !!       A = Σₘ ( yₘ + yₘ₊₁ ) · ( xₘ − xₘ₊₁ ),     m = 1…4
        !!
        !! with point 5 = point 1 to close the contour:
        !!
        !!       1 → 2 → 3 → 4 → 5 (=1)
        !!
        !! This counter-clockwise orientation yields a **positive physical area**.
        !!
        !! @note
        !! Requires the user-supplied routine
        !! `map_to_physical_coordinates(x, y, xi, eta)`.
        !!
        use mpidcl,only:decomp_info
        use MOD_solver_parameters,only:num_ghost

        class(mesh_fields), intent(inout) :: mesh
        type(decomp_info), intent(in) :: decomposition

        integer :: xi_start, eta_start, xi_end, eta_end
        integer :: i,j,m
        double precision :: face_dx,face_dy,vec_norm,area_polygon

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        !$OMP PARALLEL DO PRIVATE(m,i,face_dx,face_dy,vec_norm,area_polygon)
        do j=eta_start-num_ghost,eta_end+num_ghost
            do i=xi_start-num_ghost,xi_end+num_ghost
                call map_to_physical_coordinates(mesh%physical_space%x(i,j),mesh%physical_space%y(i,j),mesh%computational_space%xi(i,j),mesh%computational_space%eta(i,j))

                do m=1,5
                    call map_to_physical_coordinates(mesh%physical_space%x_corner(m,i,j),mesh%physical_space%y_corner(m,i,j),&
                        mesh%computational_space%xi_corner(m,i,j),mesh%computational_space%eta_corner(m,i,j))
                end do

                ! Attention: face normal vectors are pointing OUTWARDS!!

                ! TODO: add the normal vector and length ratio for the forward face!
                ! xi-sweep normal vectors
                face_dx = (mesh%physical_space%x_corner(4,i,j)-mesh%physical_space%x_corner(1,i,j))
                face_dy = (mesh%physical_space%y_corner(4,i,j)-mesh%physical_space%y_corner(1,i,j))

                vec_norm = norm2([face_dx,face_dy])
                mesh%quadrilateral_mapping%normal_vector_xi(:,1,i,j) = -face_dy/vec_norm ! x-component of normal vector
                mesh%quadrilateral_mapping%normal_vector_xi(:,2,i,j) = face_dx/vec_norm ! y-component of normal vector
                mesh%quadrilateral_mapping%length_ratio_xi(:,i,j) = vec_norm/mesh%computational_space%d_eta ! length ratio

                ! eta-sweep normal vectors
                face_dx = (mesh%physical_space%x_corner(2,i,j)-mesh%physical_space%x_corner(1,i,j))
                face_dy = (mesh%physical_space%y_corner(2,i,j)-mesh%physical_space%y_corner(1,i,j))

                vec_norm = norm2([face_dx,face_dy])
                mesh%quadrilateral_mapping%normal_vector_eta(:,1,i,j) = -face_dy/vec_norm ! x-component of normal vector
                mesh%quadrilateral_mapping%normal_vector_eta(:,2,i,j) = face_dx/vec_norm ! y-component of normal vector
                mesh%quadrilateral_mapping%length_ratio_eta(:,i,j) = vec_norm/mesh%computational_space%d_xi ! length ratio

                ! calculate the area of the cell (assume as a polygon with 4 edges -> calcualte with trapezoidal formula)
                area_polygon = 0.d0

                do m=1,4
                    area_polygon = area_polygon + 0.5*(mesh%physical_space%y_corner(m,i,j)+mesh%physical_space%y_corner(m+1,i,j))&
                        *(mesh%physical_space%x_corner(m,i,j)-mesh%physical_space%x_corner(m+1,i,j))
                end do

                ! capacity = physical area / computational area
                mesh%quadrilateral_mapping%capacity(i,j) = area_polygon / (mesh%computational_space%d_xi*mesh%computational_space%d_eta)
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine

end submodule
