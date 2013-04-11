module kinds
    implicit none
    private
    public :: sp,dp

    integer, parameter :: sp = selected_real_kind(p=4,r=30)
    integer, parameter :: dp = selected_real_kind(p=15,r=30)
end module kinds
