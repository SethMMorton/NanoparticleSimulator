/*******************************************************************
 * **********   quasi - Ellipses: n-layers
 *                      Theory:   approximate
 *                      Results:  efficiency factors
 * Axes 2 and 3 must be linked so that axis 1 moves freely and
 * axes 2 and 2 move together.  This allows us to use analytical
 * formulas for the geometrical factors.
 *
 * Seth M. Morton
 *******************************************************************/

#include "constants.h"
#include <math.h>
#include <complex.h>

int quasi (int nlayers,
           complex dielec[],
           float mdie,
           float rel_rad[][3],
           float rad[],
           float size_param,
           float *extinct,
           float *scat,
           float *absorb)
{

    /*********************
     * Perform some checks
     *********************/

    /* Too many layers for quasistatic approximatio */
    if (!(nlayers < 2)) return 1;

    /* Axes 2 and 3 must be identical */
    if (rad[1] != rad[2]) return 2;

    for (int i = 0; i < nlayers; i++) {
        /* Axes 2 and 3 must be identical */
        if (rel_rad[i][1] != rel_rad[i][2]) return 2;
    }

    *extinct = 0.0;
    *scat    = 0.0;
    *absorb  = 0.0;

    Complex(KINDR) :: die(3)
    Complex(KINDR) :: num(3), den(3)   ! Numerator and denominator

    Real(KINDR) :: tmp(nlayers,3)
    Real(KINDR) :: rel_vol(nlayers) ! Relative volume of each layer
    Real(KINDR) :: gf(3,nlayers)    ! Geometrical factors
    Real(KINDR) :: radii(3)         ! Absolute radii
    Real(KINDR) :: e                ! Eccentricity
    Real(KINDR) :: g                ! A function of eccentricity

    /*******************************
     *Calculate the relative volumes
     *******************************/

    tmp[0][:] = rel_rad[0][:];
    rel_vol[0] = PRODUCT(tmp[0][:]);
    for (int i = 1; i < nlayers; i++) {
        tmp[i][:] = SUM(rel_rad[0:i][:], DIM=1);
        rel_vol[i] = PRODUCT(tmp[i][:]) - PRODUCT(tmp[i-1][:]);
    }

    /***********************************************************
     * Calculate the goemetrical factors for each axis and layer
     ***********************************************************/
   
    for {int ilayer = 0; ilayer < nlayers; i++) {

!     Calculate the absolute radii
      radii = rel_rad(ilayer,:) * rad

!     Determine if this is a prolate or oblate spheroid or a sphere
      if (ABS(radii(1) - radii(2)) < 1E-3_KINDR) then
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

    }

!  -------------------------------
!  Determine the system dielectric
!  -------------------------------

    complex die[3];
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
    die = die / 3.0;

!  ------------------------
!  Calculate the properties
!  ------------------------

   absorb  = FOUR * size_param * AIMAG(SUM(die)) / THREE
   scat    = ( EIGHT / THREE ) * size_param**4 * ABS(SUM(die) / THREE)**2
   extinct = absorb + scat

    return 0;

}
