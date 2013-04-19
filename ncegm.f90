module ncegm
    use kinds, only: dp
    use iso_fortran_env, only: error_unit
    use ndifferential, only: lindiff
    use interpolation, only: interpolate_ordered
    use newton, only: root
    implicit none
    private
    public :: ncegm_setup,ncegm_model,ncegm_solve
    public :: ncegm_getPolicy_c,ncegm_getPolicy_aprime,ncegm_getPolicy_d,ncegm_getValueFunction
    public :: ncegm_getConditionalPolicy_c,ncegm_getConditionalPolicy_aprime,ncegm_getConditionalValueFunction
    save

    ! *******************************************************************************************************
    ! ** Types and interfaces
    ! *******************************************************************************************************
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
        integer                                                      :: a_grid_length
        procedure(ReturnFunction), pointer, nopass                   :: F=>null(),dF=>null(),d2F=>null()    ! REQUIRED: return function F and first- and second derivatives with respect to c Fc and Fcc
        procedure(ReturnFunctionMarginalInverse), pointer, nopass    :: dF_inv=>null()                      ! OPTIONAL: inverse of marginal return function Fc'
        procedure(BudgetConstraint), pointer, nopass                 :: Gamma=>null()                       ! REQUIRED: budget constraint
        procedure(BudgetConstraint), pointer, nopass                 :: dGamma=>null()                      ! OPTIONAL: derivative of budget constraint with respect to a (if present: slight increase in computation speed and accuracy
        procedure(StateTransition), pointer, nopass                  :: Psi=>null()                         ! REQUIRED if s is used: Transition function from s_t to s{t+1}
        real(dp), dimension(:,:,:), allocatable                      :: V_initial                           ! REQUIRED: initial guess of the value function
        real(dp), dimension(:), allocatable                          :: a_grid,d_grid,&                     ! REQUIRED: grids for A, D
                                                                        s_grid,z_grid                       ! REQUIRED if s or z, respectively, is used.
        real(dp)                                                     :: beta=0                              ! REQUIRED: discount factor Beta in interval (0,1)
        real(dp), dimension(:,:), allocatable                        :: z_transition                        ! REQUIRED if z is used: transition matrix for random variable z
        logical                                                      :: state_independent_foc = .FALSE.     ! OPTIONAL: is the first-order condition independent of the state variable? If set to .TRUE., then some optimizations are performed.
    end type ncegm_model

    ! *******************************************************************************************************
    ! ** Module constants
    ! *******************************************************************************************************

    integer, parameter                        :: MAX_ITER_DEFAULT = 10000               ! Default value for maximum amount of iterations
    real(dp), parameter                       :: EPSILON_DEFAULT = 0.00001_dp           ! Default value for convergence criterion epsilon
    integer, parameter                        :: STATUS_PRINT_INTERVAL_DEFAULT = 25     ! Print a status message every STATUS_PRINT_INTERVAL iteration. This is the default value.
    integer, parameter                        :: NO_STATE = -10


    ! *******************************************************************************************************
    ! ** Private model variables
    ! *******************************************************************************************************
    type(ncegm_model)                         :: model
    logical                                   :: is_setup = .FALSE., is_solved = .FALSE.
    real(dp), dimension(:,:,:), allocatable   :: valuef, policy_c, policy_aprime, policy_d
    real(dp), dimension(:,:,:,:), allocatable :: policy_c_con, policy_aprime_con, valuef_next_con
    real(dp), dimension(:,:,:,:), allocatable :: m_grid
    integer                                   :: glen_a,glen_s,glen_z,glen_d
    integer                                   :: max_iter,info_message_inverval
    real(dp)                                  :: epsilon

    contains

        ! *******************************************************************************************************
        ! ** Public subroutines to setup and solve the model
        ! *******************************************************************************************************

        ! **
        ! * Sets up and validates the input model. Note: this subroutine does not solve the model.
        ! * A further call to ncegm_solve is necessary to solve the model.
        ! *
        ! * Input:
        ! *     - dpmodel: the model to be solved
        ! *
        ! * Output:
        ! *     - status: a flag indicating if the model is setup properly. .FALSE. indicates that the model is invalid or the arrays could not be allocated.
        ! *
        ! **
        subroutine ncegm_setup(dpmodel, status)
            type(ncegm_model), intent(in) :: dpmodel
            logical, intent(out)          :: status

            status = .FALSE.
            ! Check if input model is valid
            if (.NOT. validate_model(dpmodel)) return

            ! Store model and create dummy grids if Z or S are not used
            model = dpmodel
            if (.NOT. allocated(model%s_grid)) then
                allocate(model%s_grid(1))
                model%s_grid=0.0_dp
                model%state_independent_foc = .TRUE.
            end if
            if (.NOT. allocated(model%z_grid)) then
                allocate(model%z_grid(1),model%z_transition(1,1))
                model%z_grid=0.0_dp
                model%z_transition=1.0_dp
            end if

            ! Compute grid sizes
            glen_a = size(model%a_grid)
            glen_s = size(model%s_grid)
            glen_z = size(model%z_grid)
            glen_d = size(model%d_grid)

            ! (Re)allocate all arrays
            if (allocated(valuef)) deallocate(valuef, policy_c, policy_aprime, policy_d, m_grid,policy_c_con, policy_aprime_con, valuef_next_con)
            allocate(valuef(glen_a,glen_s,glen_z), policy_c(glen_a,glen_s,glen_z), policy_aprime(glen_a,glen_s,glen_z), policy_d(glen_a,glen_s,glen_z), m_grid(glen_a,glen_d,glen_s,glen_z))
            allocate(policy_c_con(glen_a,glen_d,glen_s,glen_z), policy_aprime_con(glen_a,glen_d,glen_s,glen_z), valuef_next_con(glen_a,glen_d,glen_s,glen_z))

            ! Ready to solve model
            is_setup = .TRUE.
            status = .TRUE.
            is_solved = .FALSE.
        end subroutine ncegm_setup

        ! **
        ! * Starts the solution algorithm of the model set up by a previous call to ncegm_setup(). The algorithm terminates if a sufficient convergence is achieved
        ! * or if the maximum number of iterations is achieved. The input parameters allow to change parameters of the algorithm. The output parameter status can be used
        ! * to check if the value functions converged. All parameters are optional.
        ! *
        ! * Input:
        ! *     - max_iter_in (optional): The maximum number of value function iterations. If the value functions have not converged, the algorithm will stop after max_iter_in iterations.
        ! *     - epsilon_in (optional): The supremum norm of the difference of two subsequent value functions below which the convergence is deemed sufficient and the algorithm is stopped.
        ! *     - info_message_interval_in (optional): This value specifys the amount of iterations after which a status is printed to the output showing the current convergence.
        ! *
        ! * Output:
        ! *     - status (optional): a flag indicating if the value functions have converged
        ! **
        subroutine ncegm_solve(max_iter_in, epsilon_in, info_message_inverval_in, status)
            integer, intent(in), optional          :: max_iter_in,info_message_inverval_in
            real(dp), intent(in), optional         :: epsilon_in
            logical, intent(out), optional         :: status

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

            if (present(info_message_inverval_in)) then
                info_message_inverval = info_message_inverval_in
            else
                info_message_inverval = STATUS_PRINT_INTERVAL_DEFAULT
            end if

            ! Now solve the model
            call valfniteration(is_solved)
            if (present(status)) status = is_solved

        end subroutine ncegm_solve

        ! *******************************************************************************************************
        ! ** EGM core subroutines
        ! *******************************************************************************************************

        subroutine valfniteration(status)
            integer                                          :: index_a,index_d,index_s,index_z,index_s_prime,iter
            real(dp), dimension(glen_a,glen_s,glen_z)        :: valuef_next,dvaluef,dvaluef_next,exp_valuef,exp_dvaluef
            integer, dimension(glen_a,glen_s,glen_z)         :: d_index_choice
            real(dp)                                         :: sup_norm_diff
            real(dp), dimension(glen_a)                      :: m_end, c_end, valuef_end
            integer                                          :: k,d_max,i1
            real(dp)                                         :: d,s,z,c,t1,t2
            logical, intent(out)                             :: status

            ! 1. Compute derivative of the initial value function and the m-grid
            valuef = model%V_initial

            do index_z=1, glen_z
                do index_s=1, glen_s
                    dvaluef(:,index_s,index_z) = lindiff(model%a_grid, valuef(:,index_s, index_z))
                    do index_d=1,glen_d
                        m_grid(:,index_d,index_s,index_z) = model%Gamma(model%a_grid,model%d_grid(index_d),model%s_grid(index_s),model%z_grid(index_z))
                    end do
                end do
            end do

            call cpu_time(t1)

            ! 2. Begin value function iteration
            do iter=1,max_iter
                ! 2.1. Calculate the discounted expected value V~(a,s,z) for all a in A, s in S and z in Z
                do index_z=1,glen_z
                    do index_s=1,glen_s
                        do index_a=1,glen_a
                            exp_valuef(index_a,index_s,index_z) = model%beta*dot_product(model%z_transition(index_z,:),valuef(index_a,index_s,:))
                            exp_dvaluef(index_a,index_s,index_z) = model%beta*dot_product(model%z_transition(index_z,:),dvaluef(index_a,index_s,:))
                        end do
                    end do
                end do

                ! 2.2. Perform egm maximization and interpolation step, given (d,s,z) or (d,z) if s is not required for FOCs
                do index_z=1,glen_z
                    z=model%z_grid(index_z)
                    do index_d=1,glen_d
                        d=model%d_grid(index_d)
                        ! In case the FOCs are state-independent, s' can be computed by Psi without knowing the current state s
                        if (model%state_independent_foc) then
                            if (associated(model%Psi)) then
                                index_s_prime = model%Psi(NO_STATE,index_d,index_z)
                            else
                                index_s_prime = 1 ! Case: no state variable is used --> dummy index
                            end if
                            call egm_maximize(exp_valuef(:,index_s_prime,index_z),exp_dvaluef(:,index_s_prime,index_z),d,0.0_dp,z,m_end,c_end,valuef_end,k)

                        end if

                        do index_s=1,glen_s
                            s=model%s_grid(index_s)
                            ! The FOC depends on the state; therefore, it has to be evaluated for every state
                            if (.NOT.model%state_independent_foc) then
                                if (associated(model%Psi)) then
                                    index_s_prime = model%Psi(index_s,index_d,index_z)
                                else
                                    index_s_prime = 1
                                end if
                                call egm_maximize(exp_valuef(:,index_s_prime,index_z),exp_dvaluef(:,index_s_prime,index_z),d,s,z,m_end,c_end,valuef_end,k)
                            end if

                            ! 4. Interpolate the value function, for all (d,s,z)
                            call egm_interpolate(d,s,z,m_grid(:,index_d,index_s,index_z),m_end(1:k),c_end(1:k),valuef_end(1:k),valuef_next_con(:,index_d,index_s,index_z),policy_c_con(:,index_d,index_s,index_z))

                            policy_aprime_con(:,index_d,index_s,index_z) = m_grid(:,index_d,index_s,index_z) - policy_c_con(:,index_d,index_s,index_z)
                        end do
                    end do
                end do

                ! 3. compute unconditional value function and and its derivative
                do index_z=1,glen_z
                    z = model%z_grid(index_z)
                    do index_s=1,glen_s
                        s = model%s_grid(index_s)
                        do index_a=1,glen_a
                            d_max = maxloc(valuef_next_con(index_a,:,index_s,index_z),1)
                            d_index_choice(index_a,index_s,index_z) = d_max
                            d = model%d_grid(d_max)
                            valuef_next(index_a,index_s,index_z) = valuef_next_con(index_a,d_max,index_s,index_z)
                            if (associated(model%dGamma)) then
                                c = policy_c_con(index_a,d_max,index_s,index_z)
                                dvaluef_next(index_a,index_s,index_z) = dGamma_sc(model%a_grid(index_a),model%d_grid(d_max),model%s_grid(index_s),model%z_grid(index_z))*dF(c,d,s,z)
                            end if
                        end do
                        if (.NOT. associated(model%dGamma)) then
                            dvaluef_next(:,index_s,index_z) = lindiff(model%a_grid, valuef_next(:,index_s, index_z))
                        end if
                    end do
                end do

                ! 4. Calculate sup norm distance between last and new value function
                sup_norm_diff = maxval(abs(valuef_next-valuef))/maxval(abs(valuef))

                ! 5. Update value function
                valuef = valuef_next
                dvaluef = dvaluef_next

                ! 6. Check for convergence
                if (sup_norm_diff < epsilon) exit

                ! 7. In case the value functions have not converged yet, print a status
                if (mod(iter,info_message_inverval)==0) then
                    call cpu_time(t2)
                    print '("Iteration #",I5," sup norm distance = ",E14.6," time=",F10.3,"s")',iter,sup_norm_diff,(t2-t1)
                    t1 = t2
                end if
            end do

            if (iter==max_iter+1) then
                print *, "Error: Value functions did not converge within ",max_iter, "iterations."
                status = .FALSE.
            else
                print *, "Sufficient convergence achieved in final iteration", iter, "with supremum norm distance of", sup_norm_diff
                ! Compute unconditional policy functions
                do index_z=1,glen_z
                    do index_s=1,glen_s
                        do index_a=1,glen_a
                            d_max = d_index_choice(index_a,index_s,index_z)
                            policy_c(index_a,index_s,index_z) = policy_c_con(index_a,d_max,index_s,index_z)
                            policy_aprime(index_a,index_s,index_z) = m_grid(index_a,d_max,index_s,index_z) - policy_c(index_a,index_s,index_z)
                            policy_d(index_a,index_s,index_z)=model%d_grid(d_max)
                        end do
                    end do
                end do
                status = .TRUE.
            end if

        end subroutine

        subroutine egm_maximize(exp_valuef, exp_dvaluef, d, s, z, m_end, c_end, valuef_end, k)
            real(dp), dimension(glen_a), intent(in)            :: exp_valuef, exp_dvaluef
            real(dp), intent(in)                               :: d,s,z
            real(dp), dimension(glen_a), intent(out)           :: m_end, c_end, valuef_end
            integer, intent(out)                               :: k

            integer                                            :: index_a_prime, index_nc_a_min, index_nc_a_max, j, j_0, j_bar, max_loc
            integer,   dimension(glen_a)                       :: aprime_node

            real(dp)                                           :: bound, c_bar, exp_dval_down, exp_dval_up, v_bar, z_bar, z_temp, xguess, c_star
            real(dp), dimension(glen_a)                        :: val,c_temp
            logical, dimension(size(exp_valuef)-1)             :: mask

            index_nc_a_min = glen_a
            index_nc_a_max = 1
            exp_dval_up = 0.0_dp

            ! 1. Check for non-concavity and compute bounds of non-concave region

            ! Find discontinuities (upward jumps) in future marginal utility
            mask = (exp_dvaluef(2:glen_a)>exp_dvaluef(1:glen_a-1))
            if (count(mask)>0) then
                ! Non-concavity. Computes bounds of region.
                exp_dval_up = maxval(pack(exp_dvaluef(2:glen_a),mask))
                exp_dval_down = minval(pack(exp_dvaluef(1:glen_a-1),mask))
                if (count(exp_dvaluef>exp_dval_up)==0) then
                    index_nc_a_min=1
                else
                    index_nc_a_min = minloc(exp_dvaluef,1,exp_dvaluef>exp_dval_up)+1
                end if
                if (count(exp_dvaluef<exp_dval_down)==0) then
                    index_nc_a_max=glen_a
                else
                    index_nc_a_max = maxloc(exp_dvaluef,1,exp_dvaluef<exp_dval_down) - 1
                end if
            end if

            xguess = 1.d-1
            k = 0
            j_0 = 1
            j_bar = 0
            z_bar = 10.d6
            v_bar = exp_dval_up + 1.d-5

            do index_a_prime=1,glen_a
                if (index_a_prime<index_nc_a_min .OR. index_a_prime>index_nc_a_max) then
                    if (associated(model%dF_inv)) then
                        c_star = model%dF_inv(exp_dvaluef(index_a_prime),d,s,z)
                    else
                        c_star = root(dF_sc,d2F_sc,xguess,p2=d,p3=s,p4=z,offset_in = exp_dvaluef(index_a_prime), minx_in = 0.0_dp)
                    end if
                    k = k+1
                    c_end(k) = c_star
                    m_end(k) = c_star + model%a_grid(index_a_prime)
                    aprime_node(k) = index_a_prime
                    xguess = c_star
                else
                  if (k>0) then
                    c_bar =max(m_end(k)-model%a_grid(index_a_prime),1.d-5)
                    v_bar=dF(c_bar,d,s,z)
                    j_0 = aprime_node(k)
                  end if
                  bound = df(z_bar-model%a_grid(index_a_prime),d,s,z)

                  if (  (v_bar>exp_dvaluef(index_a_prime) .AND.  exp_dvaluef(index_a_prime).GE.bound) .OR. &
                      (exp_dvaluef(index_a_prime)<bound .AND. index_a_prime.GE.j_bar)  ) then !if_2

                    if (associated(model%dF_inv)) then
                        c_star = model%dF_inv(exp_dvaluef(index_a_prime),d,s,z)
                    else
                        c_star = root(dF_sc,d2F_sc,xguess,p2=d,p3=s,p4=z,offset_in = exp_dvaluef(index_a_prime), minx_in = 0.0_dp)
                    end if

                    z_temp = c_star + model%a_grid(index_a_prime)

                    do j = j_0,index_nc_a_max
                      c_temp(j) = z_temp - model%a_grid(j)
                      if (c_temp(j).LE.0 .AND. j>1) exit
                    end do
                    val(j_0:j-1) =  F(c_temp(j_0:j-1),d,s,z)+ exp_valuef(j_0:j-1)
                    max_loc = maxloc(val(j_0:j-1),1)+j_0-1
                    if (max_loc==index_a_prime) then
                      k = k+1
                      c_end(k) = c_star
                      m_end(k) = z_temp
                      aprime_node(k) = index_a_prime
                      xguess = c_star
                      z_bar = 10.d6
                    elseif (max_loc>index_a_prime) then
                      j_bar=max_loc
                      z_bar = z_temp
                    end if
                  end if
                end if
              end do


              valuef_end(1:k) = F(c_end(1:k),d,s,z) + exp_valuef(aprime_node(1:k))

              if (aprime_node(1)>1) then
                xguess = m_end(1)

                c_star = root(z_find_lowerbound,d_z_find_lowerbound,xguess,p2=d,p3=s,p4=z,p5=model%a_grid(aprime_node(1)),minx_in=model%a_grid(aprime_node(1)),offset_in=exp_valuef(aprime_node(1)) - exp_valuef(1))

                c_end = eoshift(c_end,-1)
                m_end = eoshift(m_end,-1)
                valuef_end = eoshift(valuef_end,-1)
                c_end(1) = c_star - 1.d-5
                m_end(1) = c_end(1) + model%a_grid(1)
                valuef_end(1) = F(c_end(1),d,s,z)+ exp_valuef(1)
                k = k + 1
              end if
        end subroutine

        subroutine egm_interpolate(d,s,z,m_exo,m_end, c_end,valuef_end,interpol_valuef,interpol_conditional_policy_c)
            real(dp), dimension(:), intent(in)        :: m_end, c_end,valuef_end
            real(dp), intent(in)                      :: d,s,z
            real(dp), dimension(glen_a), intent(in)   :: m_exo
            real(dp), dimension(glen_a), intent(out)  :: interpol_valuef,interpol_conditional_policy_c
            integer, dimension(glen_a)                :: map
            real(dp)                                  :: f_c1
            integer                                     :: i1

            f_c1 = F(c_end(1),d,s,z)
            ! Compute a mapping from the exogenous m-grid to the endogenous m-grid
            map = interpolate_ordered(m_end,m_exo)

            where (map==1) ! Constraint a'>=a_min is binding!
                where (m_exo < model%a_grid(1))
                    ! Since this case implies a nonpositive continuous choice c, and hence is infeasible, set c close to 0 and valuef as small as possible
                    interpol_conditional_policy_c = 1.d-7
                    interpol_valuef = -1.d7
                elsewhere
                    ! The constraint is binding and there is a feasible choice for c
                    interpol_conditional_policy_c = m_exo - model%a_grid(1)
                    ! Adjust the value computed function by the change in return. Note: the adjusted value function can only be less than or equal to V_end(1)
                    interpol_valuef = valuef_end(1) - f_c1 + min(F(interpol_conditional_policy_c,d,s,z),f_c1)
                end where
            elsewhere
                ! Calculate an interpolation of the consumption and value nodes on the exogenous grid of m
                interpol_conditional_policy_c = c_end(map-1)+(c_end(map)-c_end(map-1)) / (m_end(map)-m_end(map-1)) * (m_exo-m_end(map-1))
                interpol_valuef = valuef_end(map-1)+(valuef_end(map)-valuef_end(map-1)) / (m_end(map)-m_end(map-1)) * (m_exo-m_end(map-1))
            end where
        end subroutine

        ! Two callback functions used by egm_maximize to solve for the lower bound on m at and below which the constraint a >= a_min is binding
        real(dp) function z_find_lowerbound(m,d,s,z,a_prime)
            real(dp), intent(in)   :: m,d,s,z,a_prime
            z_find_lowerbound = F(m-model%a_grid(1),d,s,z) - F(m-a_prime,d,s,z)
        end function
        real(dp) function d_z_find_lowerbound(m,d,s,z,a_prime)
            real(dp), intent(in)   :: m,d,s,z,a_prime
            d_z_find_lowerbound = dF(m-model%a_grid(1),d,s,z) - dF(m-a_prime,d,s,z)
        end function




        ! *******************************************************************************************************
        ! ** Private helper subroutines
        ! *******************************************************************************************************
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
            elseif (.NOT. (allocated(m%a_grid) .AND. allocated(m%d_grid))) then
                write (error_unit,*) "Invalid model: The grids for A and D must be provided."
            else if (allocated(m%s_grid) .NEQV. associated(m%Psi)) then
                write (error_unit,*) "Invalid model: If the state variable s is used, s_grid and the transition function Psi must be provided."
            else if (allocated(m%z_grid) .NEQV. allocated(m%z_transition)) then
                write (error_unit,*) "Invalid model: If the stochastic state variable z is used, z_grid and the transition matrix must be provided."
            elseif (.NOT. allocated(m%V_initial)) then
                write (error_unit,*) "Invalid model: Missing initial guess for value function.s"
            elseif (m%beta<=0 .OR. m%beta>=1) then
                write (error_unit,*) "Invalid model: Discount factor Beta must be in interval (0,1)"
            else
                valid = .TRUE.
            end if
        end function


        ! *******************************************************************************************************
        ! ** Public getter functions for the module client to obtain the output of the algorithm
        ! *******************************************************************************************************
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

        function ncegm_getPolicy_aprime() result(aprime)
            real(dp), dimension(glen_a,glen_s,glen_z) :: aprime
            call requireSolved
            aprime = policy_aprime
        end function

        function ncegm_getValueFunction() result(v)
            real(dp), dimension(glen_a,glen_s,glen_z) :: v
            call requireSolved
            v = valuef
        end function

        function ncegm_getConditionalPolicy_c() result(c)
            real(dp), dimension(glen_a,glen_d,glen_s,glen_z) :: c
            call requireSolved
            c = policy_c_con
        end function

        function ncegm_getConditionalPolicy_aprime() result(aprime)
            real(dp), dimension(glen_a,glen_d,glen_s,glen_z) :: aprime
            call requireSolved
            aprime = policy_aprime_con
        end function

        function ncegm_getConditionalValueFunction() result(v)
            real(dp), dimension(glen_a,glen_d,glen_s,glen_z) :: v
            call requireSolved
            v = valuef_next_con
        end function

        ! *******************************************************************************************************
        ! ** Convenience functions to simplify access to array- and scalar-valued model functions
        ! *******************************************************************************************************

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
        function dGamma_sc(a,d,s,z)
            real(dp),               intent(in)  :: a
            real(dp), intent(in)                :: d,s,z
            real(dp)                             :: dGamma_sc
            real(dp), dimension(1)              :: a_arr,dg_arr
            a_arr = a
            dg_arr = model%dGamma(a_arr,d,s,z)
            dGamma_sc = dg_arr(1)
        end function
end module ncegm
