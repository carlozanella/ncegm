module ndifferential
    use kinds, only: dp
    implicit none
    private
    public :: lindiff

    contains

    ! Calculates the derivative of the function y = f(x) using the method of finite differences
    ! Input: x,y
    ! Output: the right differential df(x)/dx. Note: the last two values are equal
    function lindiff(x,y) result(dxy)
        real(dp),dimension(:), intent(in)        :: x,y
        real(dp),dimension(size(x))              :: dxy
        integer                                  :: dx_size

        dx_size = size(x)
        dxy = (eoshift(y,1)-y) / (eoshift(x,1)-x)
        dxy(dx_size)=dxy(dx_size-1)
    end function

end module ndifferential
