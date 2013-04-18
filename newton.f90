module newton
    use kinds, only: dp
    implicit none
    private
    public :: root

    contains

    real(dp) function root(f,df,xguess,minx_in,max_iter_in,eps_in, offset_in,p2,p3,p4,p5)
        procedure(real(dp))                           :: f,df
        real(dp), intent(in)                          :: xguess
        real(dp), intent(in), optional                :: p2,p3,p4,p5
        real(dp), intent(in), optional                :: eps_in, minx_in, offset_in
        integer, intent(in), optional                 :: max_iter_in

        integer                                       :: i, max_iter
        real(dp)                                      :: x,xn,dx,val,eps,minx,offset

        eps = 10.d-08
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
            if (abs(dx) <= 1.d-8*(1+abs(xn)) .AND. abs(val)<=1.d-8) then
                root = xn
                return
            end if
            x = xn
        end do
        print *, "Could not find root within ",max_iter," iterations."
        stop
    end function root
end module newton
