submodule (MOD_mathematical_functions) SMOD_eigenvalue_problem

contains

    module subroutine calculateEigenSolution(m,n,A,R,lambda_real,only_real_eigenvalues)
        !! Subroutine to calculate the eigenvalues & -vectors of a given matrix A with dimensions m x n
        !!
        !! ONLY calculates the right eigenvectors r that satisfy A * r_i = lambda_i * r_i;
        !! AND: ONLY returns the real part of the eigenvectors;
        !! Basically acts as a wrapper around LAPACK's dgeev routine
        integer, intent(in) :: m !! matrix dimension m of a m x n matrix
        integer, intent(in) :: n !! matrix dimension n of a m x n matrix
        double precision, intent(in) :: A(m,n) !! on input is the matrix A with dimensions m x n
        double precision, intent(out) :: R(m,n) !! on output is the matrix of the right eigenvectors with dimensions m x n
        double precision, intent(out) :: lambda_real(n) !! on output contains the real part of the eigenvectors Re(lambda)
        logical, intent(in) :: only_real_eigenvalues !! logical, to allow for a check at the end, if Im(lambda)=0

        integer :: lda !!leading dimension
        integer :: ldvl !! leading dimension for left eigenvectors
        integer :: ldvr !! leading dimension for right eigenvectors
        integer :: lwork !! size of work array
        integer :: info !! return info
        integer :: i,k !! loop counters

        double precision :: diag_lambda(m,n) !! at start asigned to be LCM, on output = diag(lambda)
        double precision :: lambda_imaginary(n) !! contains Im(lambda)
        double precision :: worktmp(1) !! temporary work array for query
        double precision, allocatable :: work(:) !! work array
        double precision :: dummyL(1,1) !! dummy variable for left eigenvectors, not calculated
        double precision :: dummyR(m,n) !! dummy variable for right eigenvectors, not calculated

        character(1) :: jobvl, jobvr
        jobvl = 'N'
        jobvr = 'V'

        !initialize eigenvalues
        lambda_real(:)=0.d0
        lambda_imaginary(:)=0.d0

        dummyL=0.d0
        dummyR=0.d0

        lda=n
        ldvr=n
        ldvl=n
        lwork=-1 !this is a f*** b***!! documentation in Lapack is painful..
        diag_lambda=A !create another storage for LCM, later overwritten by diag(lambda)

        ! =====================================================
        ! STEP 1: calculate the eigenvalues lambda, the matrix of right eigenvectors R
        ! =====================================================

        !query for the optimal size of work
        allocate(work(max(1,lwork)))
        worktmp(:)=0.d0
        call dgeev(jobvl,jobvr,n,diag_lambda,lda,lambda_real,lambda_imaginary,&
            dummyL,1,dummyR,ldvr,worktmp,lwork,info)

        !dynamically change the allocation for work
        lwork=int(worktmp(1))
        deallocate(work)
        allocate(work(max(1,lwork)))
        work(:)=0.d0
        !now get R and Lambda Matrices
        call dgeev(jobvl,jobvr,n,diag_lambda,lda,lambda_real,lambda_imaginary,&
            dummyL,1,R,ldvr,work,lwork,info)

        ! =====================================================
        ! STEP 2: check if one of the eigenvectors is complex (possibly violates conditions)
        ! =====================================================
        !CHECK FOR COMPLEX EIGENVALUES!!!!
        if(only_real_eigenvalues.eqv. .True.)then
            do i=1,lda
                if(lambda_imaginary(i)>ud_precision)then
                    write(*,*)"EXCEPTION OCCURED IN EIGENVALUE CALCULATION"
                    write(*,*)"-> AT LEAST ONE EIGENVALUE IS COMPLEX"
                    do k=1,lda
                        write(*,*)A(k,:)
                    end do
                    stop "ABORT"
                end if
            end do
        end if

        ! =====================================================
        ! CLEAN AND RETURN
        ! =====================================================
        deallocate(work)
    end subroutine

end submodule
