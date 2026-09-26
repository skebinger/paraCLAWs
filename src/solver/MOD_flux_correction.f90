! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_flux_correction
    !! Module containing the subroutines to perform higher order flux correction

    implicit none

contains

    subroutine limiter(decomposition,work,ixy)
        !! Limits the jump in $$ \Delta Q_{i-1/2} $$ by applying a limiter function $$ \Phi(\theta) $$ to the wave carrying the jump;<br>
        !!
        !! available flux-limiter functions: (FVMHP LeVeque equation 6.39b);<br>
        !! MINMOD: $$ \Phi(\theta) = minmod(1,\theta) $$ (FVMHP LeVeque equation 6.39b);<br>
        !! SUPERBEE: $$ \Phi(\theta) = max(0,min(1,2\theta),min(2,\theta)) $$;<br>
        !! VAN LEER: $$ \Phi(\theta) = \frac{\theta+|\theta|}{1+|\theta|} $$;<br>
        !! MC: $$ \Phi(\theta) = max(0,min((1+\theta)/2, 2, 2\theta)) $$;<br>
        !! with the slope $$ \theta = \Delta Q_{I-1/2}/\Delta Q_{i-1/2} $$, where $$ I = i-1 $$ for left going waves and $$ I = i+1 $$ for right going waves (FVMHP LeVeque equation 6.35);<br>
        !! $$ \theta $$ is calculated based on the 2-norm of the current cell i and its neighbouring cell i+1 or i-1;<br>
        !! The ratio $$ \theta_{i-1/2} $$ represents the smootheness of the data in the vicinity of i-1/2
        use mpidcl,only:decomp_info
        use MOD_work_arrays,only:work_arrays
        use MOD_solver_parameters,only:use_limiters,num_waves,num_equations,limiter_method
        use MOD_constants,only:ud_precision
        type(decomp_info), intent(in) :: decomposition
        type(work_arrays), intent(inout) :: work
        integer, intent(in) :: ixy

        !local variables:
        integer :: i,j,mw,me
        integer :: xi_start,xi_end,eta_start,eta_end
        integer :: i_right,j_right
        double precision :: wnorm2 !! 2-norm of wave
        double precision :: wnorm2_upstream
        double precision :: w_lim(num_equations)
        double precision :: r !! ratio of successive gradients
        double precision :: limiter_function !! value of the flux-limiter function

        call decomposition%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        !execute the limiter
        if(use_limiters.eqv..True.) then !skip limiter unless it is activated -> wave1d then remains unaltered

            ! precompute the right side norm at the boundary -> percompute wnorm2r(mw,xi_start-1,j) or wnorm2r(mw,i,eta_start-1)
            ! this value is chosen as the left norm beginning at i=xi_start or j=eta_start; removes need to calculate two norms per iteration
            ! helps unifying the loop for calculating the waves later on, since it needs the boundary cell for choosing the upwind cell!
            if(ixy==1)then
                i=xi_start-1
                do j=eta_start,eta_end
                    do mw=1,num_waves
                        do me=1,num_equations
                            work%wnorm2r(mw,i,j)=work%wnorm2r(mw,i,j) + work%waves(me,mw,i,j)*work%waves(me,mw,i+1,j)
                        end do
                    end do
                end do
            else
                j=eta_start-1
                do i=xi_start,xi_end
                    do mw=1,num_waves
                        do me=1,num_equations
                            work%wnorm2r(mw,i,j)=work%wnorm2r(mw,i,j) + work%waves(me,mw,i,j)*work%waves(me,mw,i,j+1)
                        end do
                    end do
                end do
            end if

            do j=eta_start,eta_end
                do i=xi_start,xi_end
                    wnorm2=0.d0
                    !loop over waves
                    do mw=1,num_waves
                        !reset to zero
                        wnorm2=0.d0
                        !assign the previous interations result for right 2-norm to left side
                        work%wnorm2l(mw,i,j)=merge(work%wnorm2r(mw,i-1,j),work%wnorm2r(mw,i,j-1),ixy==1)
                        work%wnorm2r(mw,i,j)=0.d0

                        !loop over the number of equations to calculate the 2-norm
                        !NOTE: since later, the ratio of two norms is calculated => squareroot not neccessary!
                        i_right = i + merge(1,0,ixy==1)
                        j_right = j + merge(0,1,ixy==1)
                        do me=1,num_equations
                            wnorm2=wnorm2 + work%waves(me,mw,i,j)**2
                            ! recalculate the 2-norm of the right side
                            work%wnorm2r(mw,i,j)=work%wnorm2r(mw,i,j) + work%waves(me,mw,i,j)*work%waves(me,mw,i_right,j_right)
                        end do

                        !now calculate the ratio "r" of the successive gradients
                        !here this equals to the ratio of the current wave strength
                        !and the upwind wave strength
                        wnorm2_upstream = merge(work%wnorm2l(mw,i,j),work%wnorm2r(mw,i,j),work%lambda(mw,i,j) > 0.d0)
                        r=wnorm2_upstream/wnorm2

                        select case(limiter_method)
                          case(1)!MINMOD
                            limiter_function = max(0.d0, min(1.d0,r))
                          case(2)!SUPERBEE
                            limiter_function = max(0.d0, min(1.d0, 2.d0*r), min(2.d0,r))
                          case(3)!VAN LEER
                            limiter_function = (r+abs(r)) / (1.d0+abs(r))
                          case(4)!MC
                            limiter_function = max(0.d0, min((1.d0+r)/2.d0 , 2.d0, 2.d0*r))
                          case default
                            error stop "invalid limiter choice => ABORTING"
                        end select

                        w_lim(:) = limiter_function*work%waves(:,mw,i,j)
                        !now apply the limiter to the wave in order for it to be used in the second order correction scheme
                        !choice of using merge: dont have to use cycle statement but still prevent floating point error in multiplication
                        !otherwise limiter_function would need to be 1.d0 (maybe issue with floating point arithmatic in cases where the error by multiplication
                        !gets amplified!)
                        !here: original wave is selected by merge instead of the newly computed wave w_lim
                        work%waves(:,mw,i,j)=merge(work%waves(:,mw,i,j),w_lim(:),abs(wnorm2)<=ud_precision)
                    end do
                end do
            end do
        end if
    end subroutine

    subroutine second_order_flux_correction(info,t_ctrl,work,mesh,ixy)
        !! Correction flux calculated based on capacity form differencing for high resolution methods;<br>
        !!
        !! taken from: FVMHP LeVeque, equation 6.68;<br>
        !! $$ \tilde{F}_{i-1/2}=\frac{1}{2} \sum_{p=1}^{N~waves}\left( 1-\Delta t/\Delta x / \kappa_{i-1/2} |\lambda_{i-1/2}^p| \right) |\lambda_{i-1/2}^p| \tilde{W}_{i-1/2}^p $$;<br>
        !! with: $$ \kappa_{i-1/2}=1/2 \left( \kappa_{i-1}+\kappa_i \right) $$ and<br>
        !! the modified wave $$ \tilde{W}_{i-1/2}^p = \Phi W_{i-1/2}^p $$ as the result of applying the limiter function $$ \Phi(\theta) $$
        use mpidcl,only:decomp_info
        use MOD_solver_parameters
        use MOD_time_control,only:time_control
        use MOD_work_arrays,only:work_arrays
        use MOD_domain,only:mesh_fields
        type(decomp_info), intent(in) :: info
        type(time_control), intent(in) :: t_ctrl
        type(work_arrays), intent(inout) :: work
        type(mesh_fields), intent(in) :: mesh
        integer, intent(in) :: ixy

        integer :: mw,me
        integer :: i,j
        integer :: xi_start, xi_end, eta_start, eta_end
        double precision :: delta, dt_over_delta

        call info%get_local_block_bounds(xi_start,xi_end,eta_start,eta_end)

        delta = merge(mesh%computational_space%d_xi, mesh%computational_space%d_eta, ixy==1)

        dt_over_delta = t_ctrl%get_dt()/delta

        if(ixy==1)then
            !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,me,mw)
            do j=eta_start,eta_end
                do i=xi_start-1,xi_end+1
                    do mw=1,num_waves
                        do me=1,num_equations
                            work%f(me,i,j)=work%f(me,i,j) + 0.5d0*abs(work%lambda(mw,i,j)) * &
                                ( 1.d0 - abs(work%lambda(mw,i,j)) * &
                                dt_over_delta/(0.5d0*(mesh%quadrilateral_mapping%capacity(i-1,j)+&
                                mesh%quadrilateral_mapping%capacity(i,j))) ) * work%waves(me,mw,i,j)
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        else
            !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,me,mw)
            do j=eta_start-1,eta_end+1
                do i=xi_start,xi_end
                    do mw=1,num_waves
                        do me=1,num_equations
                            work%f(me,i,j)=work%f(me,i,j) + 0.5d0*abs(work%lambda(mw,i,j)) * &
                                ( 1.d0 - abs(work%lambda(mw,i,j)) * &
                                dt_over_delta/(0.5d0*(mesh%quadrilateral_mapping%capacity(i,j-1)+&
                                mesh%quadrilateral_mapping%capacity(i,j))) ) * work%waves(me,mw,i,j)
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if
    end subroutine

end module
