!**
!* This program, like the demo_ modules, is not an integral component of this software package.
!* It only serves as a demonstration of how this package can be used.
!* This program simply runs the binary labour choice demo.
!**
program ncegm_demo
    use demo_fella11, only:start_demo_fella2011
    use demo_binary_labour_choice, only: start_demo_blc
    implicit none

    ! Uncomment the following line to solve Fella (2011)'s model
    !call start_demo_fella2011

    call start_demo_blc

end program

