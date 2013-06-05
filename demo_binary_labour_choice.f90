!**
!* This module solves a non-smooth optimization problems involving a binary labour choice using the non-concave EGM algorithm provided by ncegm.
!**
module demo_binary_labour_choice
    use kinds, only: dp
    use ncegm, only: ncegm_setup, ncegm_solve, ncegm_model, ncegm_getPolicy_c, ncegm_getPolicy_d, ncegm_getValueFunction, ncegm_getPolicy_aprime
    use grids, only: build_grid
    implicit none
    private
    public :: start_demo_blc
    save

    ! *******************************************************************************************************
    ! ** Model parameters
    ! *******************************************************************************************************
    real(dp), parameter                         :: tau    = 3_dp
    real(dp), parameter                         :: delta  = 2_dp
    real(dp), parameter                         :: roa    = 1.05_dp
    real(dp), parameter                         :: w      = 5.0_dp
    real(dp), parameter                         :: beta   = 0.95_dp


    ! *******************************************************************************************************
    ! ** Numerical parameters
    ! *******************************************************************************************************
    integer, parameter                          :: glen_a = 500
    real(dp), parameter                         :: a_min  = 0
    real(dp), parameter                         :: a_max  = 20
    integer, parameter                          :: glen_d = 2


    contains
        ! **
        ! * Sets up the household model and solves it using ncegm_solve(). Note that no state variable s (and hence no s_grid) is present in this model.
        ! **
        subroutine start_demo_blc()
            type(ncegm_model)                   :: model
            logical                             :: status
            real(dp), dimension(glen_a,1,1)     :: vfinitial,choice_d,cf,vf,ap
            integer                             :: i1

            ! Allocate required grids
            allocate(model%a_grid(glen_a),model%d_grid(glen_d))


            ! Specifies the grid for A and D
            model%a_grid = build_grid(glen_a,a_min,a_max,0) ! Uniform grid over assets a
            model%d_grid = (/ 0.0_dp , 1.0_dp /)

            ! construct a guess for the value function by assuming a consumption policy linear in available assets
            vfinitial(:,1,1) = u(0.05_dp*model%a_grid,1.0_dp,0.0_dp,0.0_dp)

            ! Specifies the functions of the model
            model%F => u
            model%dF => du
            model%d2F => d2u
            model%dF_inv => du_inv
            model%Lambda => Lambda
            model%dLambda => dLambda

            ! Specifies the discount factor Beta
            model%beta = beta

            ! Specifies an initial guess (computed above) for the value function
            model%V_initial = vfinitial

            ! Now initialize the module
            call ncegm_setup(model,status)
            if (status) then
                ! If the model is valid, the model is solved
                call ncegm_solve()
            else
                stop
            end if

            ! The computed value function and the policy functions are retrieved
            vf = ncegm_getValueFunction()
            cf = ncegm_getPolicy_c()
            ap = ncegm_getPolicy_aprime()
            choice_d = ncegm_getPolicy_d()

            ! The result is printed
            print *, "Value function V(a) & policy functions:"
            print *, ""
            print "(a12,a27,a27,a27,a27)", "a", "V(a)", "c*(a)", "a'*(a)", "d*(a)"
            print "(a12,a27,a27,a27,a27)", "----","----","----","----","----"

            do i1=1,glen_a
                print *, model%a_grid(i1), vf(i1,1,1), cf(i1,1,1), ap(i1,1,1), choice_d(i1,1,1)
            end do

        end subroutine

        ! *******************************************************************************************************
        ! ** Model functions
        ! *******************************************************************************************************

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

        function du_inv(mret,d,s,z)
            real(dp), intent(in)                :: mret,d,s,z
            real(dp)                            :: du_inv

            du_inv = (1/mret)**(1/tau)
        end function du_inv

        function d2u(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: d2u

            d2u = -tau*c**(-1-tau)

        end function d2u

        function Lambda(a,d,s,z)
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: Lambda

            Lambda = roa*a+w*d

        end function Lambda

        function dLambda(a,d,s,z)
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: dLambda

            dLambda = roa

        end function dLambda
end module demo_binary_labour_choice
