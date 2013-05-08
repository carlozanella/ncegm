!**
!* Provides functions to produce linearly and exponentially spaced grids
!**
module grids
    use kinds, only: dp
    implicit none
    private
    public :: build_grid

    contains

        ! **
        ! * Searches for the root of the function f(x, [p2,p3,p4,p5])-offset_in.
        ! *
        ! * Input:
        ! *     - grid_length: the length of the grid to be produces
        ! *     - min: the value of the first (lowest value) point in the grid
        ! *     - max: the value of the last (highest value) point in the grid
        ! *     - expg_in (optional): specifies the order of the exponential spacing of the grid. Default is 0, which
        ! *       corresponds to a linear spacing of the grid (uniform grid). A value of 1 means an exponential grid,
        ! *       a value of 2 corresponds to a double exponential grid, and so on.
        ! *
        ! * Return value: the grid
        ! **
        function build_grid(grid_length,min,max,expg_in) result(grid)
            integer, intent(in)              :: grid_length
            real(dp), intent(in)             :: min,max
            integer, intent(in), optional    :: expg_in
            integer                          :: expg,i,e
            real(dp), dimension(grid_length) :: grid
            real(dp)                         :: v1,v2,step

            if (present(expg_in)) then
                expg = expg_in
            else
                expg = 0
            end if

            v1 = min
            v2 = max

            ! Loglinearize
            do e=1,expg
                v1 = log(v1+1)
                v2 = log(v2+1)
            end do

            step = (v2-v1) / (grid_length-1)
            do i=2,grid_length
                grid(i) = v1 + (i-1)*step
            end do

            ! Exp
            do e=1,expg
                grid = exp(grid)-1
            end do

            grid(1) = min
            grid(grid_length) = max

        end function build_grid
end module
