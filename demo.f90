module demo
    use kinds, only: dp
    use ncegm, only: ncegm_model,ncegm_setup,ncegm_solve
    use grids, only: build_grid
    use tauchen, only: tauchen_discretize
    use ndifferential, only: lindiff
    implicit none
    private
    public :: start_demo



  !*********************************************************************
    ! Fella global variables
    !*********************************************************************

      INTEGER, PARAMETER   :: len_a_grid = 400   ! Asset grid size
      !----------------------------------------------------------------
      ! 1. Economic parameter (model period: half quarter)
      !----------------------------------------------------------------

      REAL (dp), PARAMETER :: theta= .77d0! Non-durables weight
      REAL (dp), PARAMETER :: tau  = 0.d0!.2435d0! Intratemporal elast.
      REAL (dp), PARAMETER :: kappa  = .075d0! Flow of house servic.
      REAL (dp), PARAMETER :: bita= 0.93d0! Discount factor

      ! AR(1) parameters for persistent income shock yp
      ! Unconditional mean mu0 normalized to 1
      REAL (dp), PARAMETER :: rho_yp = .977d0
      REAL (dp), PARAMETER :: mu0_yp = 0.d0
      REAL (dp), PARAMETER :: sigma_yp = .024d0

      ! AR(1) parameters for temporary income shock yt
      REAL (dp), PARAMETER :: rho_yt= .0d0
      REAL (dp), PARAMETER :: mu0_yt = .0d0
      REAL (dp), PARAMETER :: sigma_yt = .063d0

      REAL (dp), PARAMETER :: min_a= 0.d0! Lower bound on assets
      REAL (dp), PARAMETER :: max_a= 25.d0 ! Upper bound on assets
      REAL (dp), PARAMETER :: min_h= 1.d-2! Lower bound on housing size
      REAL (dp), PARAMETER :: max_h = 10.0  !Upper bound on housing
      REAL (dp), PARAMETER :: downp= .2d0! Downpayment ratio
      REAL (dp), PARAMETER :: gamma_= .06d0! Percentage transaction cost


      !----------------------------------------------------------------
      ! 2. Numerical parameters
      !----------------------------------------------------------------

      INTEGER, PARAMETER :: len_yp_grid = 7  ! Persistent shock grid size
      REAL(dp), PARAMETER :: cover_yp = 3.d0  ! Cover 3 SD each side

      INTEGER, PARAMETER :: len_yt_grid = 7   ! Persistent shock grid size
      REAL(dp), PARAMETER :: cover_yt = 3.d0  ! Cover 3 SD each side

      INTEGER, PARAMETER :: len_h_grid = 7  ! Number of housing choices

      REAL (dp), PARAMETER :: toler   = 1d-5  ! Numerical tolerance
      integer :: Tsimul    = 50000 ! Number of simulations for Euler errors

      ! Income grid
      INTEGER, PARAMETER :: len_y_grid = len_yp_grid*len_yt_grid


      ! Grids
      real(dp), dimension(len_a_grid), target               :: a_grid
      real(dp), dimension(len_h_grid), target               :: sd_grid
      real(dp), dimension(len_y_grid), target               :: z_grid
      real(dp), dimension(len_y_grid,len_y_grid), target    :: z_transition


      contains

        subroutine test()

        end subroutine test
        subroutine start_demo()
            type(ncegm_model) m
            real(dp), dimension(3),target :: v,w ! TODO: delete these 2 variables (here only for debugging)
            integer, parameter :: glen = 11
            real(dp), dimension(glen) :: grid
            real(dp), dimension(len_a_grid, len_h_grid, len_y_grid), target   :: vfguess
            call test()

            ! Setup the grids
            a_grid = build_grid(len_a_grid, min_a, max_a,2)
            sd_grid = build_grid(len_h_grid, min_h, max_h,0)
            call fella_create_zgrid_and_transition(z_grid, z_transition)

            ! Compute an initial guess for the value function
            vfguess = vfinitialguess()

            ! Specify return function and derivatives
            m%F=>F
            m%dF=>dF
            m%d2F=>d2F
            ! Specify budget constraint Gamma and its derivative (optional, though)
            m%Gamma => Gamma
            m%dGamma => dGamma
            ! Specify transition function for state variable s
            m%Psi=>Psi
            ! Specify Markov matrix for the transition of the stochastic variable z
            m%z_transition=>z_transition
            ! Specify the grids for A, D, S and Z
            m%a_grid=>a_grid
            m%d_grid=>sd_grid
            m%s_grid=>sd_grid
            m%z_grid=>z_grid
            ! Specify an initial guess for the value function
            m%V_initial => vfguess
            ! Finally, specify discount factor beta
            m%beta=bita
            m%state_independent_foc=.TRUE.

            call ncegm_setup(m)
            call ncegm_solve()

        end subroutine start_demo

        ! Helper functions to set up the model

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

        function vfinitialguess() result(guess)
            real(dp), dimension(len_a_grid, len_h_grid, len_y_grid)   :: guess
            integer                                                   :: i1, i2
            real(dp), dimension(len_a_grid)                           :: c_guess

            do i1 = 1,len_y_grid
              c_guess(:) =  2.d-2*a_grid(:) + z_grid(i1)
              do i2 = 1,len_h_grid
                guess(:,i2,i1) = F(c_guess, sd_grid(i2), sd_grid(i2), z_grid(i1))
              end do
            end do
        end function


        ! Functionals for the EGM algorithm

        function F(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: F
            if (tau.NE.0.d0) then
              F = LOG(theta*c**tau+(1.d0-theta)*(kappa*d)**tau)/tau
            else
              F = theta*LOG(c)+(1.d0-theta)*LOG(kappa*d)
            end if
        end function F

        function dF(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: dF

            dF = theta*c**(tau-1.d0)/(theta*c**tau+(1.d0-theta)*(kappa*d)**tau)
        end function dF

        function d2F(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: d2F

            d2F = (-tau*((theta*c**(tau-1))/(theta*c**tau+(1-theta)*(kappa*d)**tau))**2 + (tau-1)*theta*c**(tau-2)/(theta*c**tau+(1-theta)*(kappa*d)**tau))
        end function d2F

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
            dGamma = (1.d0+6.d-2)*a
        end function dGamma

        function Psi(s_index,d_index,z_index)
            integer, intent(in)                :: s_index, d_index, z_index
            integer                            :: Psi
            Psi = d_index ! new state is equal to current period's discrete choice d
        end function

end module demo
