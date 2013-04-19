! **
! * Provides a root finder using Newton's method
! **
module newton
    use kinds, only: dp
    implicit none
    private
    public :: root

    contains

    ! **
    ! * Searches for the root of the function f(x, [p2,p3,p4,p5])-offset_in.
    ! *
    ! * Input:
    ! *     - f: function, whose root will be searched
    ! *     - df: derivative of f(x,...)
    ! *     - xguess: an initial guess for the root
    ! *     - minx_in (optional): if this parameter is present, the root is only searched for in the right neighborhood of minx_in
    ! *     - max_iter_in (optional): maximum amount of iterations in Newton's method
    ! *     - eps_in (optional): if this parameter is present, the root is approximated until the absolute difference between subsequent guesses for the root abs(dx)<=eps.
    ! *       Default value, if parameter is not present: 1.d-08
    ! *     - offset_in (optional): if this parameter is present, the root is only searched for in the right neighborhood of minx_in
    ! *     - p2, p3, p4, p5 (optional): optional parameters for the function f and its derivative df
    ! *
    ! * Return value: the x* for which f(x*,...)-offset_in=0.
    ! **
    real(dp) function root(f,df,xguess,minx_in,max_iter_in,eps_in, offset_in,p2,p3,p4,p5)
        procedure(real(dp))                           :: f,df
        real(dp), intent(in)                          :: xguess
        real(dp), intent(in), optional                :: p2,p3,p4,p5
        real(dp), intent(in), optional                :: eps_in, minx_in, offset_in
        integer, intent(in), optional                 :: max_iter_in

        integer                                       :: i, max_iter
        real(dp)                                      :: x,xn,dx,val,eps,minx,offset

        eps = 1.d-08
        minx = -10.d14
        max_iter = 400
        offset = 0.0_dp
        if (present(eps_in)) eps = eps_in
        if (present(minx_in)) minx = minx_in
        if (present(max_iter_in)) max_iter = max_iter_in
        if (present(offset_in)) offset = offset_in

        x=xguess
        do i=1, max_iter
            val = f(x,p2,p3,p4,p5)-offset
            dx = val/df(x,p2,p3,p4,p5)
            xn = x - dx
            if (xn <= minx) xn=minx+0.000000000001_dp
!            if (abs(dx) <= eps .AND. xn > minx) then
            if (abs(dx) <= eps*(1+abs(xn)) .AND. abs(val)<=eps) then
                root = xn
                return
            end if
            x = xn
        end do
        print *, "Could not find root of equation within ",max_iter," iterations."
        print *, "This usually happens when the optimal continuous choice is almost 0, often indicating a problem with the model."
        stop
    end function root
end module newton
