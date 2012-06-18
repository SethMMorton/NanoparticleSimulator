program main
   
    use constants
    implicit none
    
    real(KINDR) :: extinct, scat, abs
    complex(KINDR) :: dielec(2) = (/cmplx(1.2_KINDR, 0.5_KINDR), cmplx(0.8_KINDR, 0.9_KINDR)/)
    real(KINDR) :: rel_rad(2,3) = (/0.8_KINDR, one, 0.8_KINDR, one, 0.8_KINDR, one/)
    real(KINDR) :: rad(3) = (/onefive, ten, ten/)
    call quasi(2, dielec, ONE, rel_rad, rad, THREESECOND, extinct, scat, abs)
    write(*,"('Ext ',f,' Sca ',f,' Abs ',f)") extinct, scat, abs

end program main
