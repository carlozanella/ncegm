!**
!* Module that provides a function to interpolate points on a given grid
!**
module interpolation
    use kinds, only: dp
    implicit none
    private
    public :: interpolate_ordered

    contains
        ! **
        ! * Interpolates the exogenous points on the given grid. Note: both input arrays of points are expected to be ordered.
        ! *
        ! * Input:
        ! *     - grid: array of interpolation points, i.e. the grid on which the points are interpolated
        ! *     - points: array of points that are interpolated on the given grid
        ! *
        ! * Return value: a vector mapping each given exogenous point to the corresponding RIGHT interpolating point on the grid
        ! **
        function interpolate_ordered(grid,points) result(map)
            real(dp), dimension(:), intent(in)          :: grid,points
            integer, dimension(size(points))            :: map

            integer                                     :: i,j,i_min,i_max,len_grid,len_points
            len_grid = size(grid)
            len_points = size(points)

            do i=1,len_points
                if (points(i) > grid(1)) exit
            end do
            i_min = i
            do i=len_points,1,-1
                if (points(i) < grid(len_grid)) exit
            end do
            i_max = i

            map(1:i_min-1)=1
            map(i_max+1:len_points)=len_grid

            j=2
            do i=i_min,i_max
               do while(points(i) > grid(j))
                  j = j + 1
               end do
               map(i)=j
            end do
        end function
end module interpolation




