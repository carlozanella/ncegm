module newton
    use kinds, only: dp
    implicit none
    private
    public :: root

    interface root
        module procedure root4p,root5p
    end interface

    interface
        real(dp) function functional4p(x,p2,p3,p4)
            use kinds, only: dp
            implicit none
            real(dp), intent(in)                          :: x
            real(dp), intent(in)                          :: p2,p3,p4
        end function
        real(dp) function functional5p(x,p2,p3,p4,p5)
            use kinds, only: dp
            implicit none
            real(dp), intent(in)                          :: x
            real(dp), intent(in)                          :: p2,p3,p4,p5
        end function
    end interface

    contains
    real(dp) function root4p(f,df,xguess,minx_in,max_iter_in,eps_in,offset_in,p2,p3,p4)
        procedure(functional4p)                       :: f,df
        real(dp), intent(in)                          :: xguess
        real(dp), intent(in)                          :: p2,p3,p4
        real(dp), intent(in), optional                :: eps_in, minx_in, offset_in
        integer, intent(in), optional                 :: max_iter_in

        integer                                       :: i, max_iter
        real(dp)                                      :: x,xn,dx,eps,minx,offset

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
            dx = (f(x,p2,p3,p4)-offset)/df(x,p2,p3,p4)
            xn = x - dx
            if (abs(dx) <= eps .AND. xn > minx) then
                root4p = xn
                return
            end if
            if (xn <= minx) xn=minx+0.01_dp
            x = xn
        end do
        print *, "Could not find root within ",max_iter," iterations.x=", x,"guess=",xguess
        stop

    end function root4p
    real(dp) function root5p(f,df,xguess,minx_in,max_iter_in,eps_in, offset_in,p2,p3,p4,p5)
        procedure(functional5p)                       :: f,df
        real(dp), intent(in)                          :: xguess
        real(dp), intent(in)                          :: p2,p3,p4,p5
        real(dp), intent(in), optional                :: eps_in, minx_in, offset_in
        integer, intent(in), optional                 :: max_iter_in

        integer                                       :: i, max_iter
        real(dp)                                      :: x,xn,dx,eps,minx,offset

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
            dx = (f(x,p2,p3,p4,p5)-offset)/df(x,p2,p3,p4,p5)
            xn = x - dx
            if (abs(dx) <= eps .AND. xn > minx) then
                root5p = xn
                return
            end if
            if (xn <= minx) xn=minx+0.01_dp
            x = xn
        end do
        print *, "Could not find root within ",max_iter," iterations."
        stop

    end function root5p
end module newton
