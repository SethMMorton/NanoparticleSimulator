!--------------------------------------------------------------------
! **********   mie - Spheres: n-layers
!                    Theory: exact
!                    Results: efficiency factors
! March 1999, AI SPbU
!
! Updated to FORTRAN 90 by Lasse Jensen 
! Note that the lack of documentation is due to the original code
!--------------------------------------------------------------------
Subroutine mie (nlayers, refrac_indx, rel_rad, size_param, &
                extinct, scat, absorb,                     &
                backscat, rad_pressure, albedo, asymmetry) bind ( C, name="mief")

   use iso_c_binding

   Implicit None

   Integer(C_INT),        Intent(In)  :: nlayers
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: refrac_indx(nlayers)
   Real(C_DOUBLE),    Intent(In)  :: rel_rad(nlayers)
   Real(C_DOUBLE),    Intent(In)  :: size_param
   Real(C_DOUBLE),    Intent(Out) :: extinct
   Real(C_DOUBLE),    Intent(Out) :: scat
   Real(C_DOUBLE),    Intent(Out) :: absorb
   Real(C_DOUBLE),    Intent(Out) :: backscat
   Real(C_DOUBLE),    Intent(Out) :: rad_pressure
   Real(C_DOUBLE),    Intent(Out) :: albedo
   Real(C_DOUBLE),    Intent(Out) :: asymmetry

   Integer(C_INT),        Parameter   :: nterms = 10000

   Real(C_DOUBLE)    :: d1x(nterms) 
   Complex(C_DOUBLE_COMPLEX) :: rd3x(nterms)
   Complex(C_DOUBLE_COMPLEX) :: rcx(nterms)
   Complex(C_DOUBLE_COMPLEX) :: rbb(nterms) 
   Complex(C_DOUBLE_COMPLEX) :: rd11(nterms)
   Complex(C_DOUBLE_COMPLEX) :: rd1(nterms) 
   Complex(C_DOUBLE_COMPLEX) :: rd2(nterms) 
   Complex(C_DOUBLE_COMPLEX) :: rrbb(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX) :: rrd1(nlayers,nterms) 
   Complex(C_DOUBLE_COMPLEX) :: rrd2(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX) :: srbb(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX) :: srd1(nlayers,nterms) 
   Complex(C_DOUBLE_COMPLEX) :: srd2(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX) :: ra(nterms) 
   Complex(C_DOUBLE_COMPLEX) :: rb(nterms)

   Real(C_DOUBLE)    :: xx(nlayers)

   Real(C_DOUBLE)    :: ax
   Real(C_DOUBLE)    :: ari,ari1
   Integer(C_INT)        :: i
   Integer(C_INT)        :: num, num1, num2

   Real(C_DOUBLE), Parameter :: ZERO  = 0.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: TWO   = 2.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: THREE = 3.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: FOUR  = 4.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: EIGHT = 8.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: NINE  = 9.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: HALF = ONE / TWO
   Real(C_DOUBLE), Parameter :: THIRD = ONE / THREE
   Complex(C_DOUBLE_COMPLEX), Parameter :: I_C    = CMPLX(ZERO, ONE, C_DOUBLE_COMPLEX)
   Complex(C_DOUBLE_COMPLEX), Parameter :: ZERO_C = CMPLX(ZERO, ZERO, C_DOUBLE_COMPLEX)
   Complex(C_DOUBLE_COMPLEX), Parameter :: ONE_C  = CMPLX(ONE, ZERO, C_DOUBLE_COMPLEX)

   ax = ONE / size_param
   xx(1) = size_param * rel_rad(1)
   xx(nlayers) = size_param
   do i = 2, nlayers - 1
      xx(i) = size_param * SUM(rel_rad(1:i))
   end do

!  d1(x), rd3(x), rc(x)
   num = NM(size_param)

   CALL aax(ax, num, d1x)
   CALL cd3x(size_param, num, d1x, rd3x, rcx)

   ari = ABS(refrac_indx(1))
   do i = 2, nlayers
      ari1 = ABS(refrac_indx(i))
      if(ari1.gt.ari) ari = ari1
   end do
   num2 = NM(ari*size_param)

!  rd11(m_1*x_1)
   !if(AIMAG(refrac_indx(1))*xx(1).gt.TWENTY) then
   !   write(*,*) 'k*x > 20', AIMAG(refrac_indx(1)) * xx(1)
   !end if
   CALL aa1(refrac_indx(1)*xx(1), num2, rd11)

    do i = 2, nlayers
!
!  rd1(m_i*x_i-1), rd2(m_i*x_i-1), rbb(m_i*x_i-1), rcc(m_i*x_i-1),
!
      !if (AIMAG(refrac_indx(i))*xx(i-1).gt.TWENTY) then
      !   write(*,*) 'k*x > 20', AIMAG(refrac_indx(i))*xx(i-1)
      !end if
      CALL bcd(refrac_indx(i)*xx(i-1), num2, rd1, rd2, rbb)
      rrbb(i,1:num2) = rbb(1:num2)
      rrd1(i,1:num2) = rd1(1:num2)
      rrd2(i,1:num2) = rd2(1:num2)
!
!  rd1(m_i*x_i), rd2(m_i*x_i), rbb(m_i*x_i), rcc(m_i*x_i),
!
      !if (AIMAG(refrac_indx(i))*xx(i).gt.TWENTY) then
      !   write(*,*) 'k*x > 20', AIMAG(refrac_indx(i))*xx(i)
      !end if
      CALL bcd(refrac_indx(i)*xx(i), num2, rd1, rd2, rbb)
      srbb(i,1:num2) = rbb(1:num2)
      srd1(i,1:num2) = rd1(1:num2)
      srd2(i,1:num2) = rd2(1:num2)
   end do

   CALL ABn1(nlayers, refrac_indx, num, num1, rrbb, rrd1, rrd2, &
             srbb, srd1, srd2, rd11, rd3x, rcx, d1x, ra, rb)
   CALL QQ1(ax, num1, extinct, scat, backscat, rad_pressure, ra, rb)

   absorb = extinct - scat
   albedo = scat / extinct
   asymmetry = ( extinct - rad_pressure ) / scat

Contains

!  --------------------------------------------------------------------
!  NM-auxiliary function for AA1 & BESSEL
!     (number NM is calculated using X)
!  see: Trudy Astronom. Observ. LGU V.28,P.14,1971
!     for X>1 value of NM was raised
!  August 1989, AO LGU
!  --------------------------------------------------------------------
   Function NM(x)

      Real(C_DOUBLE), Intent(In) :: x
      Integer(C_INT)                 :: NM

      if (x < 1) then
         NM = 7.5_C_DOUBLE * x + NINE
      else if (x > 100) then
         NM = 1.0625_C_DOUBLE * x + 28.5_C_DOUBLE
      else
         NM = 1.25_C_DOUBLE * x + 15.5_C_DOUBLE
      end if

   End Function NM

End Subroutine mie

! ==============
! EXTRA ROUTINES
! ==============

!  --------------------------------------------------------------------
!  AA1-subroutine for calculations of the ratio of the derivative
!    to the function for Bessel functions of half order with
!    the complex argument: J'(N)/J(N).
!    The calculations are given by the recursive expression
!    ``from top to bottom'' beginning from N=NUM.
!    RU-array of results.
!    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
!    RI - complex refractive index.
!  August 1989, AO LGU
!  --------------------------------------------------------------------
Subroutine aa1 (rx, num, ru)

   use iso_c_binding

   Implicit None

   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rx
   Integer(C_INT),        Intent(In)  :: num
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: ru(num)

   Complex(C_DOUBLE_COMPLEX) :: s, s1
   Integer(C_INT)        :: num1
   Integer(C_INT)        :: i, j
   Integer(C_INT)        :: i1 

   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE

   s = ONE / rx
   ru(num) = ( num + ONE ) * s
   num1 = num - 1
   do j = 1, num1
      i = num - j
      i1 = i + 1
      s1 = i1 * s
      ru(i) = s1 - ONE / ( ru(i1) + s1 )
   end do

End Subroutine aa1

!  --------------------------------------------------------------------
!  AAx-subroutine for calculations of the ratio of the derivative
!    to the function for Bessel functions of half order with
!    the real argument: J'(N)/J(N).
!    The calculations are given by the recursive expression
!    ``from top to bottom'' beginning from N=NUM.
!    RU-array of results.
!    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
!  March 1999, AI SPbU
!  --------------------------------------------------------------------
Subroutine AAx (a, num, ru)

   use iso_c_binding

   Implicit None

   Real(C_DOUBLE), Intent(In)  :: a
   Integer(C_INT),     Intent(In)  :: num
   Real(C_DOUBLE), Intent(Out) :: ru(num)

   Real(C_DOUBLE) :: s1
   Integer(C_INT)     :: num1
   Integer(C_INT)     :: j
   Integer(C_INT)     :: i
   Integer(C_INT)     :: i1

   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE

   ru(num) = ( num + ONE ) * a
   num1 = num - 1
   do j = 1, num1
      i = num - j
      i1 = i + 1
      s1 = i1 * a
      ru(i) = s1 - ONE / ( ru(i1) + s1 )
   end do

End Subroutine AAx

!  --------------------------------------------------------------------
!  ABn1-subroutine for calculations of the complex coefficients
!    A(N), B(N) for n-layered spheres.
!    nlayers - number of layers
!    refrac_indx(i) - complex refractive indices for innermost layer (1),
!    layer2, ... (i = 1, nlayers)
!    The coefficients are calculated up to the number NUM1.LE.NUM,
!    for which |A(N)**2+B(N)**2|.LE.10**(-40)
!    RA-array of coefficients A(N), RB-array of coefficients B(N)
!  March 1999, AI SPbU
!  --------------------------------------------------------------------
Subroutine ABn1 (nlayers, refrac_indx, num, num1, rrbb, rrd1, rrd2, &
                 srbb, srd1, srd2, rd11, rd3x, rcx, d1x,ra,rb)

   use iso_c_binding

   Implicit None

   Integer(C_INT),        Parameter   :: nterms = 10000

   Integer(C_INT),        Intent(In)  :: nlayers
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: refrac_indx(nlayers)
   Integer(C_INT),        Intent(In)  :: num
   Integer(C_INT),        Intent(Out) :: num1
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: ra(nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: rb(nterms)
   Real(C_DOUBLE),    Intent(In)  :: d1x(nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rcx(nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rd3x(nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rrbb(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rrd1(nlayers,nterms) 
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rrd2(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: srbb(nlayers,nterms)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: srd1(nlayers,nterms) 
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: srd2(nlayers,nterms)

   Complex(C_DOUBLE_COMPLEX) :: sa(nlayers) 
   Complex(C_DOUBLE_COMPLEX) :: sha(nlayers)
   Complex(C_DOUBLE_COMPLEX) :: sb(nlayers) 
   Complex(C_DOUBLE_COMPLEX) :: shb(nlayers)

   Complex(C_DOUBLE_COMPLEX) :: rd11(nterms)

   Integer(C_INT)        :: i, j

   Real(C_DOUBLE), Parameter :: ZERO  = 0.0E0_C_DOUBLE
   Complex(C_DOUBLE_COMPLEX), Parameter :: ZERO_C = CMPLX(ZERO, ZERO, C_DOUBLE_COMPLEX)

   do i = 1, num

      sa(1)  = ZERO_C
      sha(1) = rd11(i)
      sb(1)  = ZERO_C
      shb(1) = rd11(i)
!--------

      do j = 2, nlayers
         if (ABS(refrac_indx(j)*sha(j-1)-refrac_indx(j-1)*rrd2(j,i)) == ZERO) then
            sa(j) = rrbb(j,i) * (refrac_indx(j) * sha(j-1) - refrac_indx(j-1) * rrd1(j,i)) &
                 / (refrac_indx(j) * sha(j-1) - refrac_indx(j-1) * rrd2(j,i) + 1E-30_C_DOUBLE)
         else
            sa(j) = rrbb(j,i) * (refrac_indx(j) * sha(j-1) - refrac_indx(j-1) * rrd1(j,i)) &
                 / (refrac_indx(j) * sha(j-1) - refrac_indx(j-1) * rrd2(j,i))
         end if

         if (ABS(refrac_indx(j)*shb(j-1)-refrac_indx(j-1)*rrd2(j,i)).eq.ZERO) then
            sb(j) = rrbb(j,i) * (refrac_indx(j-1) * shb(j-1) - refrac_indx(j) * rrd1(j,i)) &
                 / (refrac_indx(j-1) * shb(j-1) - refrac_indx(j) * rrd2(j,i) + 1E-30_C_DOUBLE)
         else
            sb(j) = rrbb(j,i) * (refrac_indx(j-1) * shb(j-1) - refrac_indx(j) * rrd1(j,i)) &
                 / (refrac_indx(j-1) * shb(j-1) - refrac_indx(j) * rrd2(j,i))
         end if

         if (ABS(srbb(j,i) - sa(j)).eq.ZERO) then
            sha(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sa(j)) &
                  - sa(j) * srd2(j,i) / (srbb(j,i) - sa(j) + 1E-30_C_DOUBLE)
         else
            sha(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sa(j)) &
                  - sa(j) * srd2(j,i) / (srbb(j,i) - sa(j))
         end if

         if (ABS(srbb(j,i) - sb(j)).eq.ZERO) then
            shb(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sb(j)) &
                  - sb(j) * srd2(j,i) / (srbb(j,i) - sb(j) + 1E-30_C_DOUBLE)
         else
            shb(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sb(j)) &
                  - sb(j) * srd2(j,i) / (srbb(j,i) - sb(j))
         end if

      end do
!--------
! calculations of a(n), b(n)

      ra(i) = rcx(i) * (sha(nlayers) - refrac_indx(nlayers) * d1x(i)) / &
               (sha(nlayers) - refrac_indx(nlayers) * rd3x(i))

      rb(i) = rcx(i) * (refrac_indx(nlayers) * shb(nlayers) -  d1x(i)) / &
               (refrac_indx(nlayers) * shb(nlayers) - rd3x(i))

      if(ABS(ra(i))+ABS(rb(i)).le.1E-40_C_DOUBLE) exit
   end do

   num1 = i

End Subroutine ABn1

! ` --------------------------------------------------------------------
!  BCD-subroutine for calculations of the ratios of the derivative
!     to the function for Riccati-Bessel functions of half order with
!     the complex argument: psi'(N)/psi(N) and khi'(N)/khi(N)
!     and the ratios of functions: psi(N)/khi(N).
!     The calculations are given by the recursive expression
!     ``from bottom to top'' beginning from N=0.
!     rd1, rd2, rbb, rcc-arrays of results.
!     rx - (refr. index) * (size parameter)
!  March 1999, AI SPbU
!  --------------------------------------------------------------------
Subroutine bcd (rx, num, rd1, rd2, rbb)

   use iso_c_binding

   Implicit None

   Integer(C_INT),        Parameter   :: nterms = 10000

   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rx
   Integer(C_INT),        Intent(In)  :: num
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: rd1(num)
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: rd2(num)
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: rbb(num)

   Complex(C_DOUBLE_COMPLEX) :: rd3(nterms)
   Complex(C_DOUBLE_COMPLEX) :: rcc(nterms)

   Complex(C_DOUBLE_COMPLEX) :: s1
   Complex(C_DOUBLE_COMPLEX) :: rx1
   Real(C_DOUBLE)    :: x,y 

   Complex(C_DOUBLE_COMPLEX) :: rd30
   Complex(C_DOUBLE_COMPLEX) :: rxy
   Complex(C_DOUBLE_COMPLEX) :: rc0
   Complex(C_DOUBLE_COMPLEX) :: rb0
   Complex(C_DOUBLE_COMPLEX) :: r1

   Integer(C_INT) :: i

   Real(C_DOUBLE), Parameter :: ZERO  = 0.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: TWO   = 2.0E0_C_DOUBLE
   Complex(C_DOUBLE_COMPLEX), Parameter :: I_C    = CMPLX(ZERO, ONE, C_DOUBLE_COMPLEX)

   s1 = I_C
   x = REAL(rx)
   y = AIMAG(rx)
   rx1 = ONE / rx

   CALL aa1(rx, num, rd1)

!     n = 0
   rd30 = s1
   rxy = ( COS(TWO * x) + s1 * SIN(TWO * x)) * EXP(-TWO * y)
   rc0 = -( ONE - rxy ) / ( TWO * rxy )
   !rb0 = s1 * ( ONE - rxy ) / ( TWO + rxy )
   rb0 = s1 * ( ONE - rxy ) / ( ONE + rxy )
!     n = 1
   rd3(1) = -rx1 + ONE / ( rx1 - rd30 )
   rcc(1) = rc0 * ( rx1 + rd3(1) ) / ( rx1 + rd1(1) )
   rd2(1) = ( rcc(1) * rd1(1) - rd3(1) ) / ( rcc(1) - ONE )
   rbb(1) = rb0 * ( rx1 + rd2(1) ) / ( rx1 + rd1(1) )

   do i = 2, num
      r1 = i * rx1
      rd3(i) = -r1 + ONE / ( r1 - rd3(i-1) )
      rcc(i) = rcc(i-1) * ( r1 + rd3(i) ) / ( r1 + rd1(i) )
      rd2(i) = ( rcc(i) * rd1(i) - rd3(i) ) / ( rcc(i) - ONE )
      rbb(i) = rbb(i-1) * ( r1 + rd2(i) ) / ( r1 + rd1(i) )
    end do

End Subroutine bcd

!  --------------------------------------------------------------------
!  CD3X-subroutine for calculations of the ratio of the derivative
!     to the function for Riccati-Bessel functions of half order with
!     the real argument: zeta'(N)/zeta(N)
!     and the ratio of functions: psi(N)/zeta(N).
!     The calculations are given by the recursive expression
!     ``from bottom to top'' beginning from N=0.
!     rd3x, rcx-arrays of results.
!     X - size parameter
!  March 1999, AI SPbU
!  --------------------------------------------------------------------
Subroutine cd3x (x, num, d1x, rd3x, rcx)

   use iso_c_binding

   Implicit None

   Real(C_DOUBLE),    Intent(In)  :: x
   Integer(C_INT),        Intent(In)  :: num
   Real(C_DOUBLE),    Intent(In)  :: d1x(num)
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: rd3x(num)
   Complex(C_DOUBLE_COMPLEX), Intent(Out) :: rcx(num)

   Real(C_DOUBLE)    :: ax
   Real(C_DOUBLE)    :: a1

   Complex(C_DOUBLE_COMPLEX) :: rd30
   Complex(C_DOUBLE_COMPLEX) :: rxy
   Complex(C_DOUBLE_COMPLEX) :: rc0

   Integer(C_INT)        :: i

   Real(C_DOUBLE), Parameter :: ZERO  = 0.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: TWO   = 2.0E0_C_DOUBLE
   Complex(C_DOUBLE_COMPLEX), Parameter :: I_C    = CMPLX(ZERO, ONE, C_DOUBLE_COMPLEX)

   ax = ONE / x

   rd30 = I_C
   rxy = COS(TWO * x) + I_C * SIN(TWO * x)
   rc0 = - ( ONE - rxy ) / ( TWO * rxy )
   rd3x(1) = -ax + ONE / ( ax - rd30 )
   rcx(1) = rc0 * ( ax + rd3x(1) ) / ( ax + d1x(1) )

   do i = 2, num
      a1 = i * ax
      rd3x(i) = -a1 + ONE / ( a1 - rd3x(i-1) )
      rcx(i) = rcx(i-1) * ( a1 + rd3x(i) ) / ( a1 + d1x(i) )
   end do

End Subroutine cd3x

!  --------------------------------------------------------------------
!   QQ1-subroutine for calculations of the efficiency factors for
!     extinct (QEXT), scat (QSCA), backscat (QBK)
!     and radiation pressure (QPR) for spherical particles.
!   August 1989, AO LGU
!  --------------------------------------------------------------------
Subroutine QQ1 (a, num, extinct, scat, backscat, rad_pressure, ra, rb)

   use iso_c_binding

   Implicit None

   Real(C_DOUBLE),    Intent(In)  :: a
   Integer(C_INT),        Intent(In)  :: num
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: ra(num)
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: rb(num)
   Real(C_DOUBLE),    Intent(Out) :: extinct
   Real(C_DOUBLE),    Intent(Out) :: scat
   Real(C_DOUBLE),    Intent(Out) :: backscat
   Real(C_DOUBLE),    Intent(Out) :: rad_pressure

   Real(C_DOUBLE)    :: b, c, d
   Complex(C_DOUBLE_COMPLEX) :: s, r
   Integer(C_INT)        :: n, i

   Real(C_DOUBLE), Parameter :: ZERO  = 0.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: TWO   = 2.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: HALF = ONE / TWO
   Complex(C_DOUBLE_COMPLEX), Parameter :: ZERO_C = CMPLX(ZERO, ZERO, C_DOUBLE_COMPLEX)

   b = TWO * a * a
   c = ZERO
   d = ZERO
   s = ZERO_C
   r = ZERO_C
   n = 1

   do i = 1, num-1
      n = n + 2
      r = r + ( i + HALF ) * (-1)**i * ( ra(i) - rb(i) )
      s = s + i * ( i + TWO ) / ( i + ONE ) * ( ra(i) * CONJG(ra(i+1)) &
        + rb(i) * CONJG(rb(i+1)) ) + n / i / ( i + ONE )               &
        * ( ra(i) * CONJG(rb(i)) )
      c = c + n * ( REAL(ra(i)) + REAL(rb(i)) )
      d = d + n * ( ra(i) * CONJG(ra(i)) + rb(i) * CONJG(rb(i)) )
    end do

   extinct = b * c
   scat = b * d
   backscat = TWO * b * r * CONJG(r)
   rad_pressure = extinct - TWO * b * s

End Subroutine QQ1
