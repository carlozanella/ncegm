module tauchen
    use kinds, only: dp
    use normal_cdf, only: cdf_normal
    implicit none
    private
    public :: tauchen_discretize

    contains
        subroutine tauchen_discretize(y,pi,picum,rho,mu0,sigma,n,cover_tauchen)
            ! Calculates Tauchen points for AR(1) x=(1-rho)mu0 + rho*x_(-1) + u
            ! Adapted from sub_tauchen by Jesus Fernandez-Villaverde (undated)
            ! cover_tauchen is the number of SD (each side) that one wants
            ! to cover.

            ! April 10, 2013: Adapted by Carlo Zanella to use Applied Statistics
            ! algorithm as66 for computation of the normal cumulative distribution
            ! function.

            implicit none

            integer, intent(in) :: n

            real(dp), intent(in) :: rho, mu0, sigma, cover_tauchen
            real(dp), dimension (n,n), intent(out) :: pi , picum
            real(dp), dimension (n), intent(out) :: y

            integer :: yc, tc, tcc, i1

            real(dp) :: sigy, mu, upp, low

            ! Define the discrete states of the Markov chain
            sigy = sigma/sqrt(1.0-rho**2)

            y(n) = mu0 + cover_tauchen*sigy
            y(1) = mu0 - cover_tauchen*sigy

            do yc = 2, n-1

                y(yc) = y(1)+(y(n)-y(1))*(yc-1.0)/(n-1.0)

            end do

            ! Compute the transition matrix

            do tc = 1,n
                mu = (1-rho)*mu0 + rho*y(tc)
                upp = (y(2)+y(1))/2.0
                upp = (upp-mu)/sigma
                pi(tc,1) = cdf_normal(X=upp)
                low = (y(n)+y(n-1))/2.0
                low = (low-mu)/sigma
                pi(tc,n) = 1.d0-cdf_normal(X=low)

                do tcc = 2, n-1

                    low = (y(tcc-1)+y(tcc))/2.0
                    upp = (y(tcc+1)+y(tcc))/2.0
                    low = (low-mu)/sigma
                    upp = (upp-mu)/sigma
                    pi(tc,tcc) = cdf_normal(X=upp)-cdf_normal(X=low)

                end do
            end do

            ! Compute the CDF of the transition matrix.
            picum(:,1)=pi(:,1)
            do tc = 2,n
                picum(:,tc)=picum(:,tc-1)+pi(:,tc)
            enddo

        end subroutine tauchen_discretize


end module tauchen
