! **
! * Provides a demo model to be solved by the non-concave EGM algorithm. The input model is taken as is from Fella (2011).
! **
module demo_fella11
    use kinds, only: dp
    use ncegm, only: ncegm_model,ncegm_setup,ncegm_solve
    use grids, only: build_grid
    use tauchen, only: tauchen_discretize
    use ndifferential, only: lindiff
    implicit none
    private
    save
    public :: start_demo_fella2011

    ! *******************************************************************************************************
    ! ** Model parameters (Fella, 2011)
    ! *******************************************************************************************************

    integer, parameter              :: len_a_grid = 400                     ! asset grid size
    real (dp), parameter            :: theta= .77d0                         ! non-durables weight
    real (dp), parameter            :: tau  = 0.d0                          ! .2435d0! intratemporal elast.
    real (dp), parameter            :: kappa  = .075d0                      ! flow of house servic.
    real (dp), parameter            :: bita= 0.93d0                         ! discount factor

    real (dp), parameter            :: rho_yp = .977d0
    real (dp), parameter            :: mu0_yp = 0.d0
    real (dp), parameter            :: sigma_yp = .024d0

    real (dp), parameter            :: rho_yt= .0d0
    real (dp), parameter            :: mu0_yt = .0d0
    real (dp), parameter            :: sigma_yt = .063d0
    real (dp), parameter            :: downp= .2d0                          ! downpayment ratio
    real (dp), parameter            :: gamma_= .06d0                        ! percentage transaction cost


    ! *******************************************************************************************************
    ! ** Nummerical parameters
    ! *******************************************************************************************************
    real (dp), parameter            :: min_a= 0.d0                          ! lower bound on assets
    real (dp), parameter            :: max_a= 25.d0                         ! upper bound on assets
    real (dp), parameter            :: min_h= 1.d-2                         ! lower bound on housing size
    real (dp), parameter            :: max_h = 10.0                         ! upper bound on housing
    integer, parameter              :: len_yp_grid = 7                      ! persistent shock grid size
    real(dp), parameter             :: cover_yp = 3.d0                      ! cover 3 sd each side
    integer, parameter              :: len_yt_grid = 7                      ! persistent shock grid size
    real(dp), parameter             :: cover_yt = 3.d0                      ! cover 3 sd each side
    integer, parameter              :: len_h_grid = 7                       ! number of housing choices
    integer, parameter              :: len_y_grid = len_yp_grid*len_yt_grid ! Income grid


    ! *******************************************************************************************************
    ! ** Grids
    ! *******************************************************************************************************
    real(dp), dimension(len_a_grid)                       :: a_grid
    real(dp), dimension(len_h_grid)                       :: sd_grid
    real(dp), dimension(len_y_grid)                       :: z_grid
    real(dp), dimension(len_y_grid,len_y_grid)            :: z_transition

    contains

        ! **
        ! * Sets up and solves the model from Fella (2011)
        ! **
        subroutine start_demo_fella2011()
            type(ncegm_model) :: m
            logical           :: status
            real(dp), dimension(len_a_grid, len_h_grid, len_y_grid), target   :: vfguess

            ! Setup the grids
            a_grid = build_grid(len_a_grid, min_a, max_a,2)
            sd_grid = build_grid(len_h_grid, min_h, max_h,0)
            call fella_create_zgrid_and_transition(z_grid, z_transition)

            ! Compute an initial guess for the value function
            vfguess = vfinitialguess()

            ! Allocate the grids and arrays of the model
            allocate(m%a_grid(len_a_grid),m%d_grid(len_h_grid),m%s_grid(len_h_grid),m%z_grid(len_y_grid))
            allocate(m%V_initial(len_a_grid,len_h_grid,len_y_grid))
            allocate(m%z_transition(len_y_grid,len_y_grid))

            ! Specify return function and derivatives
            m%F=>u
            m%dF=>du
            m%d2F=>d2u
            ! Specify budget constraint Gamma and, optionally, its derivative
            m%Gamma => Gamma
            m%dGamma => dGamma
            ! Specify transition function for state variable s
            m%Psi=>Psi
            ! Specify Markov matrix for the transition of the stochastic variable z
            m%z_transition=z_transition
            ! Specify the grids for A, D, S and Z (note that the grids for D and S are the same)
            m%a_grid=a_grid
            m%d_grid=sd_grid
            m%s_grid=sd_grid
            m%z_grid=z_grid
            ! Specify an initial guess for the value function
            m%V_initial=vfguess
            ! Finally, specify discount factor beta
            m%beta=bita
            m%state_independent_foc=.TRUE.

            call ncegm_setup(m,status)
            if (status) call ncegm_solve()

        end subroutine start_demo_fella2011


        ! *******************************************************************************************************
        ! ** Model functions
        ! *******************************************************************************************************

        function u(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: u
            if (tau.NE.0.d0) then
              u = LOG(theta*c**tau+(1.d0-theta)*(kappa*d)**tau)/tau
            else
              u = theta*LOG(c)+(1.d0-theta)*LOG(kappa*d)
            end if
        end function u

        function du(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: du

            du = theta*c**(tau-1.d0)/(theta*c**tau+(1.d0-theta)*(kappa*d)**tau)
        end function du

        function d2u(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: d2u

            d2u = (-tau*((theta*c**(tau-1))/(theta*c**tau+(1-theta)*(kappa*d)**tau))**2 + (tau-1)*theta*c**(tau-2)/(theta*c**tau+(1-theta)*(kappa*d)**tau))
        end function d2u

        function Gamma(a,d,s,z)
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: Gamma
            real(dp)                           :: trans_cost
            if (s.NE.d) then
                trans_cost = gamma_*d
            else
                trans_cost = 0.d0
            end if
            Gamma = (1.d0+6.d-2)*( a-(1.d0-downp)*s) + s- downp*d + z - trans_cost
        end function Gamma

        function dGamma(a,d,s,z)
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: dGamma
            real(dp)                           :: trans_cost
            dGamma = 1.d0+6.d-2
        end function dGamma

        function Psi(s_index,d_index,z_index)
            integer, intent(in)                :: s_index, d_index, z_index
            integer                            :: Psi
            Psi = d_index ! new state is equal to current period's discrete choice d
        end function


        ! *******************************************************************************************************
        ! ** Helper functions to set up the model
        ! *******************************************************************************************************

        ! **
        ! * Creates the grid for the stochastic variable z consisting of the sum of two AR(1) variables discretized using Tauchen (1986)'s method.
        ! *
        ! * Output:
        ! *     - zgrid: the grid for Z
        ! *     - transition_matrix: the Markov matrix for Z
        ! **
        subroutine fella_create_zgrid_and_transition(zgrid, transition_matrix)
            real(dp), dimension(len_y_grid), intent(out)             :: zgrid
            real(dp), dimension(len_y_grid,len_y_grid), intent(out)  :: transition_matrix


            real(dp), dimension(:), allocatable                      :: yp_grid, yt_grid
            real(dp), dimension(:,:), allocatable                    :: pi_yt, picum_yt, pi_yp, picum_yp

            integer                                                  :: i1, j1, i2, j2, index

            allocate(yp_grid(len_yp_grid),pi_yp(len_yp_grid,len_yp_grid),picum_yp(len_yp_grid,len_yp_grid))
            allocate(yt_grid(len_yt_grid),pi_yt(len_yt_grid,len_yt_grid),picum_yt(len_yt_grid,len_yt_grid))

            call tauchen_discretize(yp_grid,pi_yp,picum_yp,rho_yp,mu0_yp,sigma_yp,len_yp_grid,cover_yp)
            call tauchen_discretize(yt_grid,pi_yt,picum_yt,rho_yt,mu0_yt,sigma_yt,len_yt_grid,cover_yt)

            do i1 = 1,len_yp_grid
                do j1 = 1,len_yt_grid
                    index = len_yt_grid*(i1-1) + j1
                    zgrid(index) = yp_grid(i1) + yt_grid(j1)
                end do
            end do
            zgrid = exp(zgrid)

            ! Transition matrix for z
            do i1 = 1, size(pi_yp,1)
                do i2 = 1, size(pi_yt,1)
                    do j1 = 1, size(pi_yp,2)
                        do j2 = 1, size(pi_yt,2)
                            transition_matrix((i1-1)*size(pi_yt,1)+i2,(j1-1)*size(pi_yt,2)+j2) = pi_yp(i1,j1)*pi_yt(i2,j2)
                        end do
                    end do
                end do
            end do

        end subroutine

        ! **
        ! * Computes an initial guess for the value function by assuming the household consumes 2% of its assets plus its earned income.
        ! *
        ! * Return value: initial guess for the value function
        ! **
        function vfinitialguess() result(guess)
            real(dp), dimension(len_a_grid, len_h_grid, len_y_grid)   :: guess
            integer                                                   :: i1, i2
            real(dp), dimension(len_a_grid)                           :: c_guess

            do i1 = 1,len_y_grid
              c_guess(:) =  2.d-2*a_grid(:) + z_grid(i1)
              do i2 = 1,len_h_grid
                guess(:,i2,i1) = u(c_guess, sd_grid(i2), sd_grid(i2), z_grid(i1))
              end do
            end do
        end function

end module demo_fella11
