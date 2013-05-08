! **
! * Provides the cumulative distribution function (CDF) of the normal distribution
! **
module normal_cdf
    use kinds, only: dp
    implicit none
    private
    public :: cdf_normal

    contains

        ! **
        ! * Returns the value of the normal cumulative distribution function for a given x.
        ! *
        ! * Input:
        ! *     - x: the x, at which the CDF is evaluated
        ! *
        ! * Return value: the value of the CDF of the normal distribution at x
        ! **
        real(dp) function cdf_normal(x)
            real(dp), intent(in)   :: x
            real(dp)               ::q,pdf

            call nprob(x,cdf_normal, q, pdf)
        end function cdf_normal

        subroutine nprob(z,p,q,pdf)
            real(dp), intent(in)   :: z
            real(dp), intent(out)  :: p,q,pdf
            real(dp)               :: a0,a1,a2,a3,a4,a5,a6,a7,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11
            real(dp)               :: y, zabs
        !           Reference: Adams,A.G. Areas under the normal curve,
        !           Algorithm 39, Computer J., Vol. 12, 197-8, 1969.
        !
        !           Latest Revision - 23 January 1981
        !           Obtained from: http://lib.stat.cmu.edu/apstat/
        !           April 10, 2013: modified by Carlo Zanella:
        !               - adapted to the variable naming of this program
        !               - explicit variable declarations
        !               - specified intents of dummy arguments
        !

            data a0,a1,a2,a3,a4,a5,a6,a7/0.5d0, 0.398942280444d0, &
           0.399903438504d0, 5.75885480458d0, 29.8213557808d0, &
           2.62433121679d0, 48.6959930692d0, 5.92885724438d0/, &
           b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11/0.398942280385d0, &
           3.8052d-8, 1.00000615302d0, 3.98064794d-4, 1.98615381364d0, &
           0.151679116635d0, 5.29330324926d0, 4.8385912808d0, &
           15.1508972451d0, 0.742380924027d0, 30.789933034d0, &
           3.99019417011d0/
        !
            zabs = abs(z)
            if(zabs.gt.12.7d0) go to 20
            y = a0*z*z
            pdf = exp(-y)*b0
            if(zabs.gt.1.28d0) go to 10
        !
        !       z between -1.28 and +1.28
        !
            q = a0-zabs*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
            if(z.lt.0.d0) go to 30
            p = 1.d0-q
            return
        !
        !       zabs between 1.28 and 12.7
        !
        10   q = pdf/(zabs-b1+b2/(zabs+b3+b4/(zabs-b5+b6/(zabs+b7-b8/(zabs+b9+b10/(zabs+b11))))))
            if(z.lt.0.d0) go to 30
            p = 1.d0-q
            return
        !
        !       z far out in tail
        !
        20   q = 0.d0
            pdf = 0.d0
            if(z.lt.0.d0) go to 30
            p = 1.d0
            return
        !
        !       negative z, interchange p and q
        !
        30   p = q
            q = 1.d0-p
            return
        end subroutine nprob

end module normal_cdf
