#include <math.h>
#include <complex.h>
#include "npsolve/constants.h"
#include "npsolve/mie_helpers.h"

/***************************************************
 *  NM-auxiliary function for AA1 & BESSEL
 *     (number NM is calculated using X)
 *  see: Trudy Astronom. Observ. LGU V.28,P.14,1971
 *     for X>1 value of NM was raised
 *  August 1989, AO LGU
 **************************************************/
int nm(double x) {

      if (x < 1) {
         return (int) ( 7.5 * x + 9.0 );
      } else if (x > 100) {
         return (int) ( 1.0625 * x + 28.5 );
      } else {
         return (int) ( 1.25 * x + 15.5 );
      }

}

/****************************************************************
 * AA1-subroutine for calculations of the ratio of the derivative
 *   to the function for Bessel functions of half order with
 *   the complex argument: J'(N)/J(N).
 *   The calculations are given by the recursive expression
 *   ``from top to bottom'' beginning from N=NUM.
 *   RU-array of results.
 *   A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
 *   RI - complex refractive index.
 * August 1989, AO LGU
 ****************************************************************/
void aa1 (double complex rx, int num, double complex ru[num]) {

    double complex s = 1.0 / rx;
    int num1 = num - 1;
    ru[num1] = ( num + 1 ) * s;
    
    for (int j = 0; j < num1; j++) {
        int i = num - ( j + 1 );
        int i1 = i + 1;
        double complex s1 = i1 * s;
        ru[i-1] = s1 - 1.0 / ( ru[i1-1] + s1 );
    }

}

/****************************************************************
 * AAx-subroutine for calculations of the ratio of the derivative
 *   to the function for Bessel functions of half order with
 *   the real argument: J'(N)/J(N).
 *   The calculations are given by the recursive expression
 *   ``from top to bottom'' beginning from N=NUM.
 *   RU-array of results.
 *   A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
 * March 1999, AI SPbU
 ****************************************************************/
void aax (double a, int num, double ru[num]) {

    int num1 = num - 1;
    ru[num1] = ( num + 1 ) * a;

    for (int j = 0; j < num1; j++) {
        int i = num - ( j + 1 );
        int i1 = i + 1;
        double s1 = i1 * a;
        ru[i-1] = s1 - 1.0 / ( ru[i1-1] + s1 );
    }

}

/************************************************************************
 * ABn1-subroutine for calculations of the complex coefficients
 *   A(N), B(N) for n-layered spheres.
 *   nlayers - number of layers
 *   refrac_indx(i) - complex refractive indices for innermost layer (1),
 *   layer2, ... (i = 1, nlayers)
 *   The coefficients are calculated up to the number NUM1.LE.NUM,
 *   for which |A(N)**2+B(N)**2|.LE.10**(-40)
 *   RA-array of coefficients A(N), RB-array of coefficients B(N)
 * March 1999, AI SPbU
 ************************************************************************/
int abn1 (int nlayers,
          double complex refrac_indx[nlayers],
          int num,
          int num2,
          double complex rrbb[nlayers][num2],
          double complex rrd1[nlayers][num2],
          double complex rrd2[nlayers][num2],
          double complex srbb[nlayers][num2],
          double complex srd1[nlayers][num2],
          double complex srd2[nlayers][num2],
          double complex rd11[num2],
          double complex rd3x[num],
          double complex rcx[num],
          double d1x[num],
          double complex ra[num],
          double complex rb[num]
        )
{

    double complex sa[nlayers];
    double complex sha[nlayers];
    double complex sb[nlayers];
    double complex shb[nlayers];

    int num1 = 0;
    for (int i = 0; i < num; i++) {

        sa[0]  = CMPLX(0.0, 0.0);
        sb[0]  = CMPLX(0.0, 0.0);
        sha[0] = rd11[i];
        shb[0] = rd11[i];

        for (int j = 1; j < nlayers; j++) {

            if (cabs(refrac_indx[j]*sha[j-1]-refrac_indx[j-1]*rrd2[j][i]) == 0.0) {
                sa[j] = rrbb[j][i] * ( refrac_indx[j] * sha[j-1] - refrac_indx[j-1] * rrd1[j][i] )
                      / ( refrac_indx[j] * sha[j-1] - refrac_indx[j-1] * rrd2[j][i] + 1E-30 );
            } else {
                sa[j] = rrbb[j][i] * ( refrac_indx[j] * sha[j-1] - refrac_indx[j-1] * rrd1[j][i] )
                      / ( refrac_indx[j] * sha[j-1] - refrac_indx[j-1] * rrd2[j][i] );
            }

            if (cabs(refrac_indx[j]*shb[j-1]-refrac_indx[j-1]*rrd2[j][i]) == 0.0) {
                sb[j] = rrbb[j][i] * ( refrac_indx[j-1] * shb[j-1] - refrac_indx[j] * rrd1[j][i] )
                      / ( refrac_indx[j-1] * shb[j-1] - refrac_indx[j] * rrd2[j][i] + 1E-30 );
            } else {
                sb[j] = rrbb[j][i] * ( refrac_indx[j-1] * shb[j-1] - refrac_indx[j] * rrd1[j][i] )
                      / ( refrac_indx[j-1] * shb[j-1] - refrac_indx[j] * rrd2[j][i] );
            }

            if (cabs(srbb[j][i] - sa[j]) == 0.0) {
                sha[j] = srbb[j][i] * srd1[j][i] / ( srbb[j][i] - sa[j] )
                       - sa[j] * srd2[j][i] / ( srbb[j][i] - sa[j] + 1E-30 );
            } else {
                sha[j] = srbb[j][i] * srd1[j][i] / ( srbb[j][i] - sa[j] )
                       - sa[j] * srd2[j][i] / ( srbb[j][i] - sa[j] );
            }

            if (cabs(srbb[j][i] - sb[j]) == 0.0) {
                shb[j] = srbb[j][i] * srd1[j][i] / ( srbb[j][i] - sb[j])
                       - sb[j] * srd2[j][i] / ( srbb[j][i] - sb[j] + 1E-30 );
            } else {
                shb[j] = srbb[j][i] * srd1[j][i] / ( srbb[j][i] - sb[j] )
                       - sb[j] * srd2[j][i] / ( srbb[j][i] - sb[j] );
            }

        }

        /* calculations of a(n), b(n) */
        ra[i] = rcx[i] * ( sha[nlayers-1] - refrac_indx[nlayers-1] *  d1x[i] ) /
                         ( sha[nlayers-1] - refrac_indx[nlayers-1] * rd3x[i] );

        rb[i] = rcx[i] * ( refrac_indx[nlayers-1] * shb[nlayers-1] -  d1x[i] ) /
                         ( refrac_indx[nlayers-1] * shb[nlayers-1] - rd3x[i] );

        num1++;
        if (cabs(ra[i]) + cabs(rb[i]) < 1E-40) break;

    }

    return num1;

}

/********************************************************************
 * BCD-subroutine for calculations of the ratios of the derivative
 *    to the function for Riccati-Bessel functions of half order with
 *    the complex argument: psi'(N)/psi(N) and khi'(N)/khi(N)
 *    and the ratios of functions: psi(N)/khi(N).
 *    The calculations are given by the recursive expression
 *    ``from bottom to top'' beginning from N=0.
 *    rd1, rd2, rbb, rcc-arrays of results.
 *    rx - (refr. index) * (size parameter)
 * March 1999, AI SPbU
 ********************************************************************/
void bcd (double complex rx, int num,
          double complex rd1[num], double complex rd2[num],
          double complex rbb[num])
{

    aa1(rx, num, rd1);

    double x = creal(rx);
    double y = cimag(rx);
    double complex rx1 = 1.0 / rx;

    /* n = 0 */
    double complex rd30 = I;
    double complex rxy = ( cos(2.0 * x) + I * sin(2.0 * x)) * exp(-2.0 * y);
    double complex rc0 = -( 1.0 - rxy ) / ( 2.0 * rxy );
    double complex rb0 = I * ( 1.0 - rxy ) / ( 1.0 + rxy );
    /* n = 1 */
    double complex rd3[num];
    double complex rcc[num];
    rd3[0] = -rx1 + 1.0 / ( rx1 - rd30 );
    rcc[0] = rc0 * ( rx1 + rd3[0] ) / ( rx1 + rd1[0] );
    rd2[0] = ( rcc[0] * rd1[0] - rd3[0] ) / ( rcc[0] - 1.0 );
    rbb[0] = rb0 * ( rx1 + rd2[0] ) / ( rx1 + rd1[0] );

    for (int i = 1; i < num; i++) {
        double complex r1 = ( i + 1 ) * rx1;
        rd3[i] = -r1 + 1.0 / ( r1 - rd3[i-1] );
        rcc[i] = rcc[i-1] * ( r1 + rd3[i] ) / ( r1 + rd1[i] );
        rd2[i] = ( rcc[i] * rd1[i] - rd3[i] ) / ( rcc[i] - 1.0 );
        rbb[i] = rbb[i-1] * ( r1 + rd2[i] ) / ( r1 + rd1[i] );
    }

}

/********************************************************************
 * CD3X-subroutine for calculations of the ratio of the derivative
 *    to the function for Riccati-Bessel functions of half order with
 *    the real argument: zeta'(N)/zeta(N)
 *    and the ratio of functions: psi(N)/zeta(N).
 *    The calculations are given by the recursive expression
 *    ``from bottom to top'' beginning from N=0.
 *    rd3x, rcx-arrays of results.
 *    X - size parameter
 * March 1999, AI SPbU
 ********************************************************************/
void cd3x (double x, int num, double d1x[num],
           double complex rd3x[num], double complex rcx[num]) {

    double ax = 1.0 / x;

    double complex rd30 = I;
    double complex rxy = CMPLX(cos(2.0 * x), sin(2.0 * x));
    double complex rc0 = - ( 1.0 - rxy ) / ( 2.0 * rxy );
    rd3x[0] = -ax + 1.0 / ( ax - rd30 );
    rcx[0]  = rc0 * ( ax + rd3x[0] ) / ( ax + d1x[0] );

    for (int i = 1; i < num; i++) {
        double a1 = ( i + 1 ) * ax;
        rd3x[i] = -a1 + 1.0 / ( a1 - rd3x[i-1] );
        rcx[i] = rcx[i-1] * ( a1 + rd3x[i] ) / ( a1 + d1x[i] );
    }

}

/****************************************************************
 *  QQ1-subroutine for calculations of the efficiency factors for
 *    extinct (QEXT), scat (QSCA), backscat (QBK)
 *    and radiation pressure (QPR) for spherical particles.
 *  August 1989, AO LGU
 ****************************************************************/
void qq1 (double a, int num, int num1, double *extinct, double *scat,
          double *backscat, double *rad_pressure,
          double complex ra[num], double complex rb[num])
{

    double b = 2.0 * SQR(a);
    double c = 0.0;
    double d = 0.0;
    double complex s = CMPLX(0.0, 0.0);
    double complex r = CMPLX(0.0, 0.0);
    int n = 1;

    //printf("%d, %d\n", num, num1);
    for (int i = 0; i < num1-1; i++) {
        //printf("%f, %f, %f, %f\n", creal(ra[i]), cimag(ra[i]), creal(rb[i]), cimag(rb[i]));
        int i1 = i + 1;
        n += 2;
        r += ( i1 + 0.5 ) * pow(-1.0, i1) * ( ra[i] - rb[i] );
        s += i1 * ( i1 + 2.0 ) / ( i1 + 1.0 ) * ( ra[i] * conj(ra[i+1])
           + rb[i] * conj(rb[i+1]) ) + n / i1 / ( i1 + 1.0 )
           * ( ra[i] * conj(rb[i]) );
        c += n * ( creal(ra[i]) + creal(rb[i]) );
        d += n * ( ra[i] * conj(ra[i]) + rb[i] * conj(rb[i]) );
    }

    *extinct = b * c;
    *scat = b * d;
    *backscat = 2.0 * b * r * conj(r);
    *rad_pressure = *extinct - 2.0 * b * s;

}
