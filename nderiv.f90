! **
! * Provides functions to compute the derivative of a function numerically
! **
module nderiv
    use kinds, only: dp
    implicit none
    private
    public :: linderiv

    contains

        ! **
        ! * Calculates the derivative of the function y = f(x) using the method of finite differences.
        ! *
        ! * Input:
        ! *     - x: array representing the x-values of the function f
        ! *     - y: array representing the y-values of the function f (same size as array x)
        ! *
        ! * Return value: array for the right derivatives df(x)/dx.
        ! *
        ! * Note: the last two values are equal
        ! **
        function linderiv(x,y) result(dxy)
            real(dp),dimension(:), intent(in)        :: x,y
            real(dp),dimension(size(x))              :: dxy
            integer                                  :: dx_size

            dx_size = size(x)
            dxy = (eoshift(y,1)-y) / (eoshift(x,1)-x)
            dxy(dx_size)=dxy(dx_size-1)
        end function

end module nderiv
