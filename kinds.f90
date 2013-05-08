!**
!* Provides kinds used by all modules of this package.
!* This module can be understood as an interface to integrate this package in other projects.
!**
module kinds
    implicit none
    private
    public :: dp

    integer, parameter :: dp = selected_real_kind(p=15,r=30)
end module kinds
