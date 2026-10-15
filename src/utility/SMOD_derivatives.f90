submodule(MOD_mathematical_functions) SMOD_derivatives

contains

    module function deriv_2p(x1,x2,y1,y2) result(deriv)
        !! compute the derivative of the function y=f(x) from 2 points
        double precision :: x1, x2, y1, y2
        double precision :: deriv

        deriv = (y2-y1)/(x2-x1)
    end function

end submodule
