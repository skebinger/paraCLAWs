submodule(MOD_mathematical_functions) SMOD_matrix_operations

contains

    module subroutine matrix_vector_product(m,n,A,v,b)
        !! Subroutine to calculate the simple matrix vector product B=A*v
        integer, intent(in) :: m !! matrix dimension m of a m x n matrix
        integer, intent(in) :: n !! matrix dimension n of a m x n matrix
        double precision,intent(in) :: A(m,n) !! on input is the matrix A with dimensions m x n
        double precision,intent(in) :: v(n) !! vector v; with v element of |R^n
        double precision,intent(out) :: b(n) !! on output is the product result b=A*v, with dimension n

        !local variables
        integer :: lda !! leading matrix dimension
        double precision :: alpha, beta !alpha=1, beta=0; hard coded!

        lda=m
        alpha=1.0
        beta=0.0

        !calculate: B = A * x
        call dgemv('N', m, n, alpha, A ,lda ,v ,1 ,beta , b, 1)
    end subroutine

    module subroutine inverse_matrix(m,n,A)
        !! Subroutine to calculate the inverse of a matrix A
        !!
        !! On input: A = matrix to be inverted; On output: inverse_A = inv(A). <br>
        !! Acts as a wrapper around LAPACK's dgetrf and dgetri routines
        integer, intent(in) :: m !! matrix dimension m of a m x n matrix
        integer, intent(in) :: n !! matrix dimension n of a m x n matrix
        double precision,intent(inout) :: A(m,n) !! on input is the matrix A with dimensions m x n; on output is the inerse of A with dimensions m x n

        !local variables
        integer :: lda !! leading matrix dimension
        integer :: lwork !! size of work array
        integer :: info !! return info
        double precision, allocatable :: work(:) !! work array
        integer :: ipiv(m,n) !! pivot index matrix of LU-Factorization

        lda=m
        lwork=-1 !this is a f*** b***!! documentation in Lapack is painful..
        allocate(work(1))
        work(:)=0.d0

        ! =====================================================
        ! STEP 1: calculate the LU-factorization
        ! =====================================================
        !on output: A contains the factors L and U
        call dgetrf(m,n,A,lda,ipiv,info)

        ! =====================================================
        ! STEP 2: calculate the inv(A) based on the LU factorization
        ! =====================================================
        !first get a query for the sice of work
        call dgetri(n,A,lda,ipiv,work,lwork,info)
        lwork=int(work(1))
        if(allocated(work))then
            deallocate(work)
            allocate(work(lwork))
        else
            stop "allocation error in alpha calculation"
        end if
        !now use the pivot matrix and LU-factoriz. to invert A
        call dgetri(n,A,lda,ipiv,work,lwork,info)

    end subroutine

end submodule
