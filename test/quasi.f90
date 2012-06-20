!--------------------------------------------------------------------
! **********   quasi - Ellipses: n-layers
!                      Theory:   approximate
!                      Results:  efficiency factors
! Axes 2 and 3 must be linked so that axis 1 moves freely and 
! axes 2 and 2 move together.  This allows us to use analytical
! formulas for the geometrical factors.
!
! Seth M. Morton
!--------------------------------------------------------------------
!Subroutine quasi (nlayers, dielec, mdie, rel_rad, rad, size_param, &
Subroutine quasi (nlayers, dielec, mdie, rr, rad, size_param, &
                  extinct, scat, absorb) bind ( C, name="quasif" )

   use iso_c_binding

   Implicit None

   Integer(C_INT),        Intent(In)  :: nlayers            ! Number of layers
   Complex(C_DOUBLE_COMPLEX), Intent(In)  :: dielec(nlayers)    ! Dielectric of layers
   Real(C_DOUBLE),    Intent(In)  :: mdie               ! Dielectric of medium
   !Real(C_DOUBLE),    Intent(In)  :: rel_rad(nlayers,3) ! Relative layer radius
   Real(C_DOUBLE),    Intent(In)  :: rr(nlayers*3) ! Relative layer radius
   Real(C_DOUBLE),    Intent(In)  :: rad(3)             ! Radius of each axis
   Real(C_DOUBLE),    Intent(In)  :: size_param         ! Size parameter
   Real(C_DOUBLE),    Intent(Out) :: extinct            ! Extiction
   Real(C_DOUBLE),    Intent(Out) :: scat               ! Scattering
   Real(C_DOUBLE),    Intent(Out) :: absorb             ! Absorption

   Complex(C_DOUBLE_COMPLEX) :: die(3)
   Complex(C_DOUBLE_COMPLEX) :: num(3), den(3)   ! Numerator and denominator

   Real(C_DOUBLE) :: tmp(nlayers,3)
   Real(C_DOUBLE) :: rel_vol(nlayers) ! Relative volume of each layer
   Real(C_DOUBLE) :: gf(3,nlayers)    ! Geometrical factors
   Real(C_DOUBLE) :: radii(3)         ! Absolute radii
   Real(C_DOUBLE) :: e                ! Eccentricity
   Real(C_DOUBLE) :: g                ! A function of eccentricity

   Integer(C_INT)     :: i, ilayer

   Real(C_DOUBLE), Parameter :: ZERO  = 0.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: ONE   = 1.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: TWO   = 2.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: THREE = 3.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: FOUR  = 4.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: EIGHT = 8.0E0_C_DOUBLE
   Real(C_DOUBLE), Parameter :: THIRD = ONE / THREE
   Real(C_DOUBLE), Parameter :: PI    = 3.14159265358979323846_C_DOUBLE

   Real(C_DOUBLE) :: rel_rad(nlayers,3)
   rel_rad = RESHAPE(rr, (/nlayers, 3/))

   if (nlayers > 2) then
      write(*,*) 'Too many layers for quasistatic approximation'
      stop
   end if

   if (rad(2) /= rad(3)) then
      write(*,*) 'Axes 2 and 3 must be identical'
      stop
   end if

   do i = 1, nlayers
      if (rel_rad(i,2) /= rel_rad(i,3)) then
         write(*,*) 'Axes 2 and 3 must be identical'
         stop
      end if
   end do

   extinct = ZERO
   scat    = ZERO
   absorb  = ZERO

!  ------------------------------
!  Calculate the relative volumes
!  ------------------------------

   tmp(1,:) = rel_rad(1,:)
   rel_vol(1) = PRODUCT(tmp(1,:))
   do i = 2, nlayers
      tmp(i,:) = SUM(rel_rad(1:i,:), DIM=1)
      rel_vol(i) = PRODUCT(tmp(i,:)) - PRODUCT(tmp(i-1,:))
   end do

!  ---------------------------------------------------------
!  Calculate the goemetrical factors for each axis and layer
!  ---------------------------------------------------------
   
   do ilayer = 1, nlayers

!     Calculate the absolute radii
      radii = rel_rad(ilayer,:) * rad

!     Determine if this is a prolate or oblate spheroid or a sphere
      if (ABS(radii(1) - radii(2)) < 1E-3_C_DOUBLE) then
!        Sphere
         gf(:,ilayer) = THIRD
      else if (radii(1) > radii(2)) then
!        Prolate (cigar)
!        Determine long axis
         e = ONE - ( radii(2) / radii(1) )
         gf(1,ilayer) = ( ( ONE - e ) / e )                           &
                      * ( -ONE + ( ONE / ( 2 * SQRT(e) ) )            &
                         * LOG(( ONE + SQRT(e) ) / ( ONE - SQRT(e) )) &
                        )
!        Set short axes.  Total must be 1, and the short are equal
         gf(2:3,ilayer) = ( ONE - gf(1,ilayer) ) / TWO
      else
!        Oblate (pancake)
!        Determine long axes
         e = ONE - ( radii(1) / radii(2) )
         g = SQRT(( ONE - e ) / e)
         gf(2:3,ilayer) = ( g / ( TWO * e ) ) * ( ( PI / TWO ) - ATAN(g) ) &
                        - ( g**2 / TWO )
!        Set short axis.  Total must be 1, and the long are equal
         gf(1,ilayer) = ONE - SUM(gf(2:3,ilayer))
      end if

   end do

!  -------------------------------
!  Determine the system dielectric
!  -------------------------------

!  One layer
   if (nlayers == 1) then
      die = ( dielec(1) - mdie ) / ( mdie + gf(:,1) * ( dielec(1) - mdie ) )
!  Two layers
   else if (nlayers == 2) then
!     Numerator
      num = ( dielec(1) - dielec(2) ) * ( gf(:,1) - rel_vol(1) * gf(:,2) )
      num = ( dielec(2) - mdie ) * ( dielec(2) + num )
      num = num + rel_vol(1) * dielec(2) * ( dielec(1) - dielec(2) )
!     Denominator
      den = ( dielec(1) - dielec(2) ) * ( gf(:,1) - rel_vol(1) * gf(:,2) )
      den = ( dielec(2) + den  ) * ( mdie + ( dielec(2) - mdie ) * gf(:,2) )
      den = den + rel_vol(1) * gf(:,2) * dielec(2) * ( dielec(1) - dielec(2) )
!     The whole thing
      die = num / den
   end if
!  Factor of three omitted above for ease
   die = die / THREE

!  ------------------------
!  Calculate the properties
!  ------------------------

   absorb  = FOUR * size_param * AIMAG(SUM(die)) / THREE
   scat    = ( EIGHT / THREE ) * size_param**4 * ABS(SUM(die) / THREE)**2
   extinct = absorb + scat

End Subroutine quasi
