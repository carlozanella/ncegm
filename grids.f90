module grids
    use kinds, only: dp
    implicit none
    private
    public :: build_grid

    contains
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
