ncegm
=====

What is it?
-----------
ncegm (non-concave EGM) is a software for solving non-concave dynamic
programming models with one continuous and one discrete choice variable.
The implementation is based on the algorithm by Fella (2011), which is a
generalization of the endogenous gridpoints method (EGM) by Caroll (2006)
to non-concave problems.

Usage
-----
The module `ncegm` contains all functions to solve a specific dynamic
programming model. The basic usage is the following:

1. Import the necessary functions from the ncegm module.
2. Define the input model and allocate the grids.
3. Call ncegm_setup(model) to initialize ncegm with the given model.
4. Call ncegm_solve() to solve the model.
5. Use getter functions as documented below to obtain the results.

### ncegm functions
The functions and its parameters are documented in detail in the source code itself.
The following provides an overview of all public functions of ncegm.

* ncegm_setup: sets up the module for a specific input model
* ncegm_solve: solves the model previously set up by ncegm_setup

The following functions can be called after the model is solved to obtain the results.

* ncegm_getPolicy_c: the policy function for c for a given state (a,s,z)
* ncegm_getPolicy_aprime: the policy function for a' for a given state (a,s,z)
* ncegm_getPolicy_d: the policy function for d for a given state (a,s,z)
* ncegm_getValueFunction: the value function

These functions have the same interpretation as above, except that they return the optimal
policies and their values conditional on a choice d.

* ncegm_getConditionalPolicy_c
* ncegm_getConditionalPolicy_aprime
* ncegm_getConditionalValueFunction

### Type ncegm_model and interfaces
An input model is specified by an instance of the derived type `ncegm_model`.

    type ncegm_model
        procedure(ReturnFunction), pointer, nopass                   :: F=>null(),dF=>null()
        procedure(ReturnFunction), pointer, nopass                   :: d2F=>null()
        procedure(ReturnFunctionMarginalInverse), pointer, nopass    :: dF_inv=>null()
        procedure(BudgetConstraint), pointer, nopass                 :: Lambda=>null()
        procedure(BudgetConstraint), pointer, nopass                 :: dLambda=>null()
        procedure(StateTransition), pointer, nopass                  :: Psi=>null()
        real(dp), dimension(:,:,:), allocatable                      :: V_initial
        real(dp), dimension(:), allocatable                          :: a_grid,d_grid,&
                                                                        s_grid,z_grid
        real(dp)                                                     :: beta=0
        real(dp), dimension(:,:), allocatable                        :: z_transition
        logical                                                      :: state_independent_foc = .FALSE.
    end type ncegm_model

Note that there are 4 procedure interfaces, which are used by the model.

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

Not all fields are required to specify a complete model. The elements of
the model are documented in the following sections:

### Return function and budget constraint
The following functions are required:

* F(c,d,s,z): is the return function
* dF(c,d,s,z): is the first derivative of the return function
* Lambda(a,s,z,d): is the budget constraint function

The following function is optional:

* dLambda(a,s,z,d): derivative of the budget constraint function with respect to a
           This is used to apply the envelope theorem to compute V'(a,s,z).
           It is recommended to provide this function.

### Marginal inverse of the return function
These functions are used to find the numerical root of the first-order condition. Only one
of them have to be provided.

* d2F(c,d,s,z): the second-order derivative of the return function
* dF_inv(mret,d,s,z): the marginal inverse of the return function

### Transition function / matrix
These functions are used if discrete state variables are present in the model.

* Psi(s,d,z): deterministic transition function of the model. If this function is present,
       then it is assumed that s is used by the model. Therefore, a grid for S is required.
* z_transition: transition matrix for stochastic variable z. If this matrix is present,
       then it is assumed that z is used by the model. Therefore, a grid for Z is requried.

### Grids
To use ncegm, grids for the state and discrete decision variables have to be supplied. Grids
are assumed to be ordered. Note that the grids must be allocated before their use.

The following grids are always required:

* a_grid: grid for the continuous state variable a. Note that a_grid(1) is a_min.
* d_grid: grid for the discrete choice variable d.

The following grids are required if s, z or both are used, respectively.

* s_grid: grid for the deterministic state variable s. Required iff Psi is present.
* z_grid: grid for the stochastic state variable z. Required iff z_transition is present.

### Initial guess for V
* V_initial: this is the initial guess for the value function. Make sure its numerical derivative
       is not zero. Note that the dimension of V are A*S*Z where A, S and Z are the
       sizes of the respective grids.

### Further parameters
* state_independent_foc (boolean flag): if s is used by the model but the first-order conditions do not depend on
       s, a significant speed up can be achieved by setting this flag to true. The first-order
       conditions do not depend on s if neither dF(c,d,s,z) nor Psi(s,d,z) depend on s. An example
       where this is the case is in Fella (2011)'s model, where s'=d.

### Further options
For further options, check out the documentation in the source code. ncegm_solve provides some additional configuration
options.

Example usage
-------------
A simple example is provided in demo_binary_labour_choice.f90. The following snippet illustrates
how the model is defined. Note that no s_grid and z_grid are defined because neither s nor z is
used as a state variable in this simple model.

    use ncegm, only: ncegm_setup, ncegm_solve, ncegm_model, ncegm_getPolicy_c, ncegm_getPolicy_d, ncegm_getValueFunction, ncegm_getPolicy_aprime

    (...)

        type(ncegm_model)                   :: model

        (...)

        ! Allocate required grids
        allocate(model%a_grid(glen_a),model%d_grid(glen_d))

        (...)

        model%a_grid = build_grid(glen_a,a_min,a_max,0) ! Uniform grid over assets a
        model%d_grid = (/ 0.0_dp , 1.0_dp /)

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
        call ncegm_setup(model)
        (...)
        call ncegm_solve()

        ! The computed value function and the policy functions are retrieved
        vf = ncegm_getValueFunction()
        cf = ncegm_getPolicy_c()
        ap = ncegm_getPolicy_aprime()
        choice_d = ncegm_getPolicy_d()
        (...)

This model is already implemented and is solved by the produced binary as a demonstration (see main.f90).

Building ncegm
--------------
ncegm was developed and tested with gfortran but also builds fine with g95.
The target standard is Fortran 2003 (although many of its features are used
because they are not yet implemented in gfortran).

Warnings regarding unused variables in demo_binary_labour_choice.f90 and
demo_fella11.f90 come from the fact that the model functions do not use
all model variables - either because they do not exist or because they do not
depend on it for other reasons.

To compile with g95 the provided Makefile has to be edited. Note that some flags need to be translated when compiling
with g95.

Further notes
-------------
When trying to solve a model with a huge state space, the program might run out of stack (Segmentation fault.). This can be solved in one of two ways:

* either disable optimizations so that arrays are no longer allocated on the stack (a performance penalty has
  to be expected), or
* increase the stack size: ulimit -s [new-stack-size] with bash on Linux.

Todo
----
* Arbitrary dimensionality of the state and decision space
* Provide an option to solve finite horizon problems.
  Although this requires that value functions are saved after
  every iteration.

Licensing
---------
See LICENSE.
