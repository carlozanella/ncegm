module ncegm
    use kinds, only: dp
    use iso_fortran_env, only: error_unit
    use ndifferential, only: lindiff
    implicit none
    private
    public :: ncegm_setup,ncegm_model,ncegm_solve
    save

    ! Type and interface declarations
    interface
        function ReturnFunction(c,d,s,z) result(r)
            use kinds, only: dp
            implicit none
            real(dp), dimension(:), intent(in) :: c
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(c))       :: r
        end function
        function ReturnFunctionMarginalInverse(mret,d,s,z) result(c)
            use kinds, only: dp
            implicit none
            real(dp), intent(in)               :: mret,d,s,z
            real(dp)                           :: c
        end function
        function BudgetConstraint(a,d,s,z) result(m)
            use kinds, only: dp
            implicit none
            real(dp), dimension(:), intent(in) :: a
            real(dp), intent(in)               :: d,s,z
            real(dp), dimension(size(a))       :: m
        end function
        function StateTransition(s_index,d_index,z_index) result(s_prime_index)
            use kinds, only: dp
            implicit none
            integer, intent(in)                :: s_index, d_index, z_index
            integer                            :: s_prime_index
        end function
    end interface

    interface F
        module procedure F_sc,F_arr
    end interface F
    interface dF
        module procedure dF_sc,dF_arr
    end interface dF
    interface d2F
        module procedure d2F_sc,d2F_arr
    end interface d2F


    type ncegm_model
        procedure(ReturnFunction), pointer, nopass                   :: F=>null(),dF=>null(),d2F=>null()    ! REQUIRED: return function F and first- and second derivatives with respect to c Fc and Fcc
        procedure(ReturnFunctionMarginalInverse), pointer, nopass    :: dF_inv=>null()                      ! OPTIONAL: inverse of marginal return function Fc'
        procedure(BudgetConstraint), pointer, nopass                 :: dGamma=>null()                      ! OPTIONAL: derivative of budget constraint with respect to a (if present: slight increase in computation speed and accuracy
        procedure(BudgetConstraint), pointer, nopass                 :: Gamma=>null()                       ! REQUIRED: budget constraint
        procedure(StateTransition), pointer, nopass                  :: Psi=>null()                         ! REQUIRED: Transition function from s_t to s{t+1}
        real(dp), dimension(:,:,:), pointer                          :: V_initial=>null()                   ! REQUIRED: initial guess of the value function
        real(dp), dimension(:), pointer                              :: a_grid=>null(),s_grid=>null(),&     ! REQUIRED: grids for A, D, S and Z
                                                                        z_grid=>null(),d_grid=>null()
        real(dp)                                                     :: beta=0                              ! REQUIRED: discount factor Beta in interval (0,1)
        real(dp), dimension(:,:), pointer                            :: z_transition=>null()                ! REQUIRED: transition matrix for random variable z
    end type ncegm_model

    ! Module constants
    integer, parameter                        :: MAX_ITER_DEFAULT = 10000       ! Default value for maximum amount of iterations
    real(dp), parameter                       :: EPSILON_DEFAULT = 0.00001_dp   ! Default value for convergence criterion epsilon


    ! Private module variables
    type(ncegm_model)                         :: model
    logical                                   :: is_setup = .FALSE., is_solved = .FALSE.
    real(dp), dimension(:,:,:), allocatable   :: valuef, policy_c, policy_d
    real(dp), dimension(:,:,:,:), allocatable :: m_grid
    integer                                   :: glen_a,glen_s,glen_z,glen_d
    integer                                   :: max_iter
    real(dp)                                  :: epsilon

    contains
        ! Public module functions to set up the model
        subroutine ncegm_setup(dpmodel)
            type(ncegm_model), intent(in) :: dpmodel
            ! TODO: check if model is valid
            if (.NOT. validate_model(dpmodel)) return




            ! Store model and grid sizes
            model = dpmodel
            glen_a = size(model%a_grid)
            glen_s = size(model%s_grid)
            glen_z = size(model%z_grid)
            glen_d = size(model%d_grid)
            ! (Re)allocate all arrays
            if (allocated(valuef)) deallocate(valuef, policy_c, policy_d, m_grid)
            allocate(valuef(glen_a,glen_s,glen_z), policy_c(glen_a,glen_s,glen_z), policy_d(glen_a,glen_s,glen_z), m_grid(glen_a,glen_d,glen_s,glen_z))
            is_setup = .TRUE.
            is_solved = .FALSE.
        end subroutine ncegm_setup

        subroutine ncegm_solve(max_iter_in, epsilon_in)
            integer, intent(in), optional          :: max_iter_in
            real(dp), intent(in), optional         :: epsilon_in

            if (.NOT.is_setup) return

            ! Prepare configurations of the algorithm

            if (present(max_iter_in)) then
                max_iter = max_iter_in
            else
                max_iter = MAX_ITER_DEFAULT
            end if

            if (present(epsilon_in)) then
                epsilon = epsilon_in
            else
                epsilon = EPSILON_DEFAULT
            end if

            ! Now solve the model
            call valfniteration
            is_solved = .TRUE.
        end subroutine ncegm_solve

        ! Public module getter methods
        function ncegm_getPolicy_d() result(d)
            real(dp), dimension(glen_a,glen_s,glen_z) :: d
            call requireSolved
            d = policy_d
        end function

        function ncegm_getPolicy_c() result(c)
            real(dp), dimension(glen_a,glen_s,glen_z) :: c
            call requireSolved
            c = policy_c
        end function

        function ncegm_getValueFunction() result(v)
            real(dp), dimension(glen_a,glen_s,glen_z) :: v
            call requireSolved
            v = valuef
        end function

        ! EGM algorithm functions

        subroutine valfniteration()
            integer                                   :: index_d,index_s,index_z,iter
            real(dp), dimension(glen_a,glen_s,glen_z) :: valuef_next,dvaluef,exp_valuef,exp_dvaluef
            real(dp)                                  :: sup_norm_diff

            ! 1. Compute derivative of initial value function and m-grid
            valuef = model%V_initial
            do index_z=1, glen_z
                do index_s=1, glen_s
                    dvaluef(:,index_z,index_s) = lindiff(model%a_grid, valuef(:,index_z, index_s))
                    do index_d=1,glen_d
                        m_grid(:,index_d,index_s,index_z) = model%Gamma(model%a_grid,model%d_grid(index_d),model%s_grid(index_s),model%z_grid(index_z))
                    end do
                end do
            end do

            ! 1b) Compute m-grid
            do index_z=1,glen_z
                do index_s=1,glen_s
                    do index_d=1,glen_d
                        m_grid(:,index_d,index_s,index_z) = model%Gamma(model%a_grid,model%d_grid(index_d),model%s_grid(index_s),model%z_grid(index_z))
                    end do
                end do
            end do

            do iter=1,max_iter
                ! 2. Solve for discounted expected value V~
                !--
                ! 3. perform egm step
                !--
                ! 4. interpolate valuefn
                !--
                ! 5. compute unconditional value function
                !--
                valuef = 1
                valuef_next = 1
                ! 6. Check for convergence
                sup_norm_diff = MAXVAL(ABS(valuef_next-valuef))/MAXVAL(ABS(valuef))
                if (sup_norm_diff < epsilon) exit
                ! 7. Update value function
                valuef = valuef_next
            end do

            if (iter==max_iter) then
                print *, "Value functions did not converge."
            else
                print *, ""
            end if

        end subroutine

        subroutine egm_maximize(exp_valuef, exp_dvaluef, d, s, z)
            real(dp), dimension(glen_a), intent(in)            :: exp_valuef, exp_dvaluef
            real(dp), intent(in)                               :: d, s, z

            ! TODO: implement:
                !-> also output: consumption choice and m-values!!



        end subroutine

        subroutine egm_interpolate()

        end subroutine







        ! Some convenience functions
        function F_sc(c,d,s,z)
            real(dp), intent(in)   :: c,d,s,z
            real(dp)               :: F_sc
            real(dp), dimension(1) :: c_arr, ret
            c_arr=c
            ret = F(c_arr,d,s,z)
            F_sc = ret(1)
        end function F_sc
        function F_arr(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: F_arr
            F_arr=model%F(c,d,s,z)
        end function F_arr
        function dF_sc(c,d,s,z)
            real(dp), intent(in)   :: c,d,s,z
            real(dp)               :: dF_sc
            real(dp), dimension(1) :: c_arr, ret
            c_arr=c
            ret = dF(c_arr,d,s,z)
            dF_sc = ret(1)
        end function dF_sc
        function dF_arr(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: dF_arr
            dF_arr=model%dF(c,d,s,z)
        end function dF_arr
        function d2F_sc(c,d,s,z)
            real(dp), intent(in)   :: c,d,s,z
            real(dp)               :: d2F_sc
            real(dp), dimension(1) :: c_arr, ret
            c_arr=c
            ret = d2F(c_arr,d,s,z)
            d2F_sc = ret(1)
        end function d2F_sc
        function d2F_arr(c,d,s,z)
            real(dp), dimension(:), intent(in)  :: c
            real(dp), intent(in)                :: d,s,z
            real(dp), dimension(size(c))        :: d2F_arr
            d2F_arr=model%d2F(c,d,s,z)
        end function d2F_arr

        ! Helper subroutines
        subroutine requireSolved()
            if (.NOT.is_solved) then
                print *, "Error: The model must first be solved before obtaining policy and value functions."
                stop
            end if
        end subroutine

        logical function validate_model(m) result(valid)
            type(ncegm_model), intent(in) :: m
            valid = .FALSE.
            if (.NOT. (associated(m%F) .AND. associated(m%dF) .AND. associated(m%d2F))) then
                write (error_unit,*) "Invalid model: The return function and first and second derivatives with respect to c must be provided."
            elseif (.NOT. associated(m%Gamma)) then
                write (error_unit,*) "Invalid model: The budget constraint function Gamma must be provided"
            elseif (.NOT. associated(m%Psi)) then
                write (error_unit,*) "Invalid model: The state variable transition function Psi must be provided"
            elseif (.NOT. (associated(m%a_grid) .AND. associated(m%s_grid) .AND. associated(m%z_grid) .AND. associated(m%d_grid))) then
                write (error_unit,*) "Invalid model: The grids for A, S,Z and D must be provided."
            elseif (.NOT. associated(m%z_transition)) then
                write (error_unit,*) "Invalid model: Missing transition matrix for random variable z"
             elseif (.NOT. associated(m%V_initial)) then
                write (error_unit,*) "Invalid model: Missing initial guess for value function.s"
            elseif (m%beta<=0 .OR. m%beta>=1) then
                write (error_unit,*) "Invalid model: Discount factor Beta must be in interval (0,1)"
            else
                valid = .TRUE.
            end if
        end function
end module ncegm
