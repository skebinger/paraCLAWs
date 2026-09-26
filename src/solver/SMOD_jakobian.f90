! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

submodule (MOD_flux_calculation) SMOD_jakobian
    !!Submodule containing jakobian related calculations
    implicit none

contains

    module function getJakobi(ixy,i,j,mesh,solution,work,decomposition)result(jakobi)
        !! Returns the jakobian matrix for the state vars(i,j); w/ aux array
        use mpi_f08
        use MOD_solver_parameters,only:num_equations,psi=>grav
        use MOD_user_subroutines
        use MOD_user_data
        use MOD_index
        use MOD_domain,only:solution_fields,mesh_fields
        use MOD_mathematical_functions, only: integrate
        use MOD_work_arrays,only:work_arrays
        use mpidcl,only:decomp_info
        integer, intent(in) :: ixy !! direction of sweep
        integer, intent(in) :: i,j
        type(mesh_fields), intent(in) :: mesh
        type(solution_fields), intent(in) :: solution
        type(work_arrays), intent(in) :: work
        type(decomp_info), intent(in) :: decomposition
        double precision :: jakobi(num_equations,num_equations)

        double precision :: top(6)

        integer:: k,l
        !double precision :: psi = 0.d0
        integer :: rank,ierr

        double precision :: h,qu,qv
        double precision :: Ui,Vi,Uw,Vw
        double precision :: deltaU,deltaV
        double precision :: beta,F,G,dF,dG
        double precision :: dA1dh,dA1dqu,dA1dqv
        double precision :: dA2dh,dA2dqu,dA2dqv
        double precision :: dA3dh,dA3dqu,dA3dqv
        double precision :: dB1dh,dB1dqu,dB1dqv
        double precision :: dB2dh,dB2dqu,dB2dqv
        double precision :: dB3dh,dB3dqu,dB3dqv
        integer :: region_ident
        !==========================================================
        ! START OF USER INPUT
        !==========================================================

        !if(ixy==1)then
        !    jakobi(1,:) = [0.d0, 1.d0, 0.d0]
        !    jakobi(2,:) = [-vars(2,i,j)**2/vars(1,i,j)**2 + psi*vars(1,i,j), 2.d0*vars(2,i,j)/vars(1,i,j), 0.d0]
        !    jakobi(3,:) = [-vars(2,i,j)*vars(3,i,j)/vars(1,i,j)**2, vars(3,i,j)/vars(1,i,j), vars(2,i,j)/vars(1,i,j)]
        !else
        !    !stop
        !    jakobi(1,:) = [0.d0, 0.d0, 1.d0]
        !    jakobi(2,:) = [-vars(2,i,j)*vars(3,i,j)/vars(1,i,j)**2, vars(3,i,j)/vars(1,i,j), vars(2,i,j)/vars(1,i,j)]
        !    jakobi(3,:) = [-vars(3,i,j)**2/vars(1,i,j)**2 + psi*vars(1,i,j), 0.d0, 2.d0*vars(3,i,j)/vars(1,i,j)]
        !end if

        jakobi=0.d0

        do k=1,num_equations
            do l=1,num_equations
                if(isnan(jakobi(k,l)))then
                    write(*,*)"isnan in jakobi"
                    write(*,*) "cell: ", i,j
                    write(*,*) "direction: ", ixy
                    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
                    print *, "rank ", rank
                    write(*,*) solution%vars(:,i,j)
                    write(*,*) jakobi(1,:),1
                    write(*,*) jakobi(2,:),2
                    write(*,*) jakobi(3,:),3
                    error stop
                end if
            end do
        end do
    end function

end submodule
