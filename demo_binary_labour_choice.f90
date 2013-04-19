module demo_binary_labour_choice
    use kinds, only: dp
    use ncegm, only: ncegm_setup, ncegm_solve, ncegm_model, ncegm_getPolicy_c, ncegm_getPolicy_d, ncegm_getValueFunction, ncegm_getPolicy_aprime
    use grids, only: build_grid
    private
    public :: start_demo_blc
    save

    ! Model parameters
    real(dp), parameter                         :: tau    = 0.2435_dp
    real(dp), parameter                         :: delta  = 1.75_dp
    real(dp), parameter                         :: roa    = 1.05_dp
    real(dp), parameter                         :: w      = 3.0_dp
    real(dp), parameter                         :: beta   = 0.93_dp

    ! Nummerical parameters
    integer, parameter                          :: glen_a = 400
    real(dp), parameter                         :: a_min  = 0
    real(dp), parameter                         :: a_max  = 100
    integer, parameter                          :: glen_d = 2


    contains
        subroutine start_demo_blc()
            type(ncegm_model)                   :: model
            real(dp), dimension(glen_a,1,1)     :: vfinitial,choice_d,cf,vf,ap
            integer                             :: i1

            ! Allocate required grids
            allocate(model%a_grid(glen_a),model%d_grid(glen_d))


            model%a_grid = build_grid(glen_a,a_min,a_max,0) ! Uniform grid over assets a
            model%d_grid = (/ 0.0_dp , 1.0_dp /)


            vfinitial(:,1,1) = u(0.1_dp*model%a_grid,1.0_dp,0.0_dp,0.0_dp)


            model%F => u
            model%dF => du
            model%d2F => d2u
            model%Gamma => Gamma
            model%dGamma => dGamma
            model%beta = beta
            model%V_initial = vfinitial

            call ncegm_setup(model)
            call ncegm_solve()
            vf = ncegm_getValueFunction()
            cf = ncegm_getPolicy_c()
            ap = ncegm_getPolicy_aprime()
            choice_d = ncegm_getPolicy_d()
            print *, "Value function V(a) & policy functions:"
            print *, ""
            print "(a12,a27,a27,a27,a27)", "a", "V(a)", "c*(a)", "a'*(a)", "d*(a)"
            print "(a12,a27,a27,a27,a27)", "----","----","----","----","----"

            do i1=1,glen_a
                print *, model%a_grid(i1), vf(i1,1,1), cf(i1,1,1), ap(i1,1,1), choice_d(i1,1,1)
            end do

        end subroutine

        function u(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: u

            u =  c**(1-tau)/(1-tau) - delta*d

        end function u

        function du(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: du

            du = c**(-tau)

        end function du

        function d2u(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: d2u

            d2u = -tau*c**(-1-tau)

        end function d2u

        function Gamma(a,d,s,z)
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: Gamma

            Gamma = roa*a+w*d

        end function Gamma

        function dGamma(a,d,s,z)
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: dGamma

            dGamma = roa

        end function dGamma
end module demo_binary_labour_choice
