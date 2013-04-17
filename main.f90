program ncegm_demo
    use demo_fella11, only:start_demo_fella2011
    use demo_binary_labour_choice, only: start_demo_blc
    implicit none

!    call start_demo_fella2011
    call start_demo_blc

end program

