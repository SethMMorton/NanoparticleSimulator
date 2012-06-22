/***************************************************************
* **********   mie - Spheres: n-layers
*                    Theory: exact
*                    Results: efficiency factors
* March 1999, AI SPbU
*
* Updated to FORTRAN 90 by Lasse Jensen 
* Converted to C from Fortran 90 by Seth M. Morton
* Note that the lack of documentation is due to the original code
*****************************************************************/

#include <math.h>
#include <complex.h>
#include "npsolve/constants.h"
#include "npsolve/arraycontractions.h"
#include "npsolve/solvers.h"

#define SQR(x) ((x) * (x))
#define CMPLX(r, i) ((r) + I*(i))

/**************************************
 * Declarations of supporting functions
 **************************************/

int nm(double x);

void aa1 (double complex rx, int num, double complex ru[num]);

void aax (double a, int num, double ru[num]);

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
        );

void bcd (double complex rx, int num,
          double complex rd1[num], double complex rd2[num],
          double complex rbb[num]);

void cd3x (double x, int num, double d1x[num],
           double complex rd3x[num], double complex rcx[num]);

void qq1 (double a, int num, int num1, double *extinct, double *scat,
          double *backscat, double *rad_pressure,
          double complex ra[num], double complex rb[num]);

/****************************
 * The main Mie theory solver 
 ****************************/

int mie (int nlayers,                           /* Number of layers */
         double complex refrac_indx[MAXLAYERS], /* Refrac. index of layers */
         double rel_rad[MAXLAYERS],             /* Relative radii of layers */
         double size_param,                     /* Size parameter */
         double *extinct,                       /* Extinction */
         double *scat,                          /* Scattering */
         double *absorb,                        /* Absorption */
         double *backscat,                      /* Backscattering */
         double *rad_pressure,                  /* Radiation pressure */
         double *albedo,                        /* Albedo */
         double *asymmetry                      /* Asymmetry */
        )
{

    /* Return code... normally 0, but 1 if outside of reasonable range */
    int retcode = 0;

    double xx[nlayers];
    double ax = 1.0 / size_param;
    xx[0] = size_param * rel_rad[0];
    xx[nlayers-1] = size_param;
    for (int i = 1; i < nlayers - 1; i++) {
        xx[i] = size_param * sumrange(nlayers, rel_rad, 0, i);
    }

    /* d1(x), rd3(x), rc(x) */
    int num = nm(size_param);

    double d1x[num];
    aax(ax, num, d1x);
    double complex rd3x[num], rcx[num];
    cd3x(size_param, num, d1x, rd3x, rcx);

    double ari = cabs(refrac_indx[0]);
    for (int i = 1; i < nlayers; i++) {
        double tmp = cabs(refrac_indx[i]);
        if (tmp > ari) ari = tmp;
    }
    int num2 = nm(ari*size_param);

    /* rd11(m_1*x_1) */
    if (cimag(refrac_indx[0]) * xx[0] > 20.0) {
        /* k*x > 20 AIMAG(refrac_indx(1)) * xx(1) */
        retcode = 1;
    }
    double complex rd11[num2];
    aa1(refrac_indx[0]*xx[0], num2, rd11);

    double complex rbb[num2], rd1[num2], rd2[num2];
    double complex rrbb[nlayers][num2], rrd1[nlayers][num2];
    double complex rrd2[nlayers][num2], srbb[nlayers][num2];
    double complex srd1[nlayers][num2], srd2[nlayers][num2];
    for (int i = 1; i < nlayers; i++) {

        /* rd1(m_i*x_i-1), rd2(m_i*x_i-1), rbb(m_i*x_i-1), rcc(m_i*x_i-1) */
        if (cimag(refrac_indx[i]) * xx[i-1] > 20.0) {
            /* k*x > 20 AIMAG(refrac_indx(i))*xx(i-1) */
            retcode = 1;
        }
        bcd(refrac_indx[i]*xx[i-1], num2, rd1, rd2, rbb);
        for (int j = 0; j < num2; j++) {
            rrbb[i][j] = rbb[j];
            rrd1[i][j] = rd1[j];
            rrd2[i][j] = rd2[j];
        }

        /* rd1(m_i*x_i), rd2(m_i*x_i), rbb(m_i*x_i), rcc(m_i*x_i) */
        if (cimag(refrac_indx[i]) * xx[i] > 20.0) {
            /* k*x > 20 AIMAG(refrac_indx(i))*xx(i) */
            retcode = 1;
        }
        bcd(refrac_indx[i]*xx[i], num2, rd1, rd2, rbb);
        for (int j = 0; j < num2; j++) {
            srbb[i][j] = rbb[j];
            srd1[i][j] = rd1[j];
            srd2[i][j] = rd2[j];
        }
    }

    double complex ra[num], rb[num];
    int num1 = abn1(nlayers, refrac_indx, num, num2, rrbb, rrd1, rrd2,
                    srbb, srd1, srd2, rd11, rd3x, rcx, d1x, ra, rb);
    qq1(ax, num, num1, extinct, scat, backscat, rad_pressure, ra, rb);

    *absorb = *extinct - *scat;
    *albedo = *scat / *extinct;
    *asymmetry = ( *extinct - *rad_pressure ) / *scat;

    return retcode;

}

/*************************************
 * Definitions of supporting functions
 *************************************/

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
