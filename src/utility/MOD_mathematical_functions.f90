! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_mathematical_functions
    !! Module containing some useful mathematical functions or methods
    use MOD_constants,only:ud_precision
    implicit none

    abstract interface
        function f_of_x(x)result(fx)
            double precision :: x
            double precision :: fx
        end function
    end interface

    abstract interface
        function f_of_x_and_g(x,g)result(fx)
            import::f_of_x
            double precision :: x
            procedure(f_of_x) :: g
            double precision :: fx
        end function
    end interface

    ! interface for derivatives
    interface
        module function deriv_2p(x1,x2,y1,y2) result(deriv)
            !! compute the derivative of the function y=f(x) from 2 points
            double precision :: x1, x2, y1, y2
            double precision :: deriv
        end function
    end interface

    ! interface for the eigenvalue calculation in jakobian
    interface
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
        end subroutine
    end interface

    ! interface for the matrix operations
    interface
        module subroutine matrix_vector_product(m,n,A,v,b)
            !! Subroutine to calculate the simple matrix vector product B=A*v
            integer, intent(in) :: m !! matrix dimension m of a m x n matrix
            integer, intent(in) :: n !! matrix dimension n of a m x n matrix
            double precision,intent(in) :: A(m,n) !! on input is the matrix A with dimensions m x n
            double precision,intent(in) :: v(n) !! vector v; with v element of |R^n
            double precision,intent(out) :: b(n) !! on output is the product result b=A*v, with dimension n
        end subroutine
        module subroutine inverse_matrix(m,n,A)
            !! Subroutine to calculate the inverse of a matrix A
            !!
            !! On input: A = matrix to be inverted; On output: inverse_A = inv(A). <br>
            !! Acts as a wrapper around LAPACK's dgetrf and dgetri routines
            integer, intent(in) :: m !! matrix dimension m of a m x n matrix
            integer, intent(in) :: n !! matrix dimension n of a m x n matrix
            double precision,intent(inout) :: A(m,n) !! on input is the matrix A with dimensions m x n; on output is the inerse of A with dimensions m x n
        end subroutine
    end interface

contains

    function exp_func(x,gamma,sigma)result(fx)
        double precision :: x
        double precision :: fx
        double precision :: gamma
        double precision :: sigma

        fx=exp(-(x-1)**gamma/sigma)
    end function

    function square_func(x,g)result(fx)
        double precision :: x
        double precision :: fx
        procedure(f_of_x) :: g

        fx=g(x)*g(x)
    end function

    function heavyside(x)result(h)
        !! Heavyside step function
        !!
        !! H(x)=0, if x<0
        !!
        !! H(x)=1, if x>=0
        double precision, intent(in) :: x !! function argument
        double precision :: h !! result H(x) of step function (1 or 0)

        if(x>=0) then
            h=1.d0
        else
            h=0.d0
        end if
    end function

    function integrate(f, g, b, lb, ub, eps, max_iter,adaptive) result(integral)
        !! Function to calculate the integral of a function f(x) or g(x) over the interval [lb, ub]
        !!
        !! If f is provided, it will be used for integration; otherwise, g will be used.
        !! f is defind as a function of x and g(x), while g is defined as a function of x only.
        !!
        !! If b is provided, it will be used as a scaling factor for the integration variable
        !!
        implicit none

        double precision, intent(in) :: lb, ub
        double precision, intent(in), optional :: eps
        integer, intent(in) :: max_iter
        logical :: adaptive

        procedure(f_of_x_and_g), optional :: f
        procedure(f_of_x) :: g
        double precision, optional, intent(in) :: b

        double precision :: integral, prev_integral
        double precision :: h, k
        integer :: n, i, iter

        if (present(b)) then
            k = b
        else
            k = 1.d0
        end if

        if (.not. adaptive) then
            ! ==================================================
            ! Fixed integration
            ! ==================================================
            n = max_iter
            h = (ub - lb) / n
            integral = 0.d0

            if (present(f)) then
                do i = 0, n-1
                    integral = integral + f((lb + i*h)*k, g) * h
                end do
            else
                do i = 0, n-1
                    integral = integral + g((lb + i*h)*k) * h
                end do
            end if

        else
            ! ==================================================
            ! Adaptive integration with relative tolerance
            ! ==================================================
            n = 50
            integral = 0.d0
            prev_integral = huge(1.d0)

            iter = 0
            do while (iter < max_iter)

                prev_integral = integral
                integral = 0.d0
                h = (ub - lb) / n

                if (present(f)) then
                    do i = 0, n-1
                        integral = integral + f((lb + i*h)*k, g) * h
                    end do
                else
                    do i = 0, n-1
                        integral = integral + g((lb + i*h)*k) * h
                    end do
                end if

                ! Relative convergence check
                if (iter > 0) then
                    if (abs(integral - prev_integral) / &
                        max(1.d0, abs(prev_integral)) < eps) then
                        !print *, "Integral relative convergence: ", abs(integral - prev_integral) / max(1.d0, abs(prev_integral))
                        !print *, "Integral value: ", integral
                        exit
                    end if
                end if

                n = int(n * 1.5)
                iter = iter + 1

            end do
        end if

    end function integrate


end module
