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

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "npsolve/constants.h"
#include "npsolve/arraycontractions.h"
#include "npsolve/solvers.h"
#include "npsolve/mie_helpers.h"

int mie (int nlayers,                         /* Number of layers */
         double complex refrac_indx[MAXLAYERS], /* Refractive index of layers */
         double rel_rad[MAXLAYERS],             /* Relative radii of layers */
         double size_param,                   /* Size parameter */
         double *extinct,                     /* Extinction */
         double *scat,                        /* Scattering */
         double *absorb,                      /* Absorption */
         double *backscat,                    /* Backscattering */
         double *rad_pressure,                /* Radiation pressure */
         double *albedo,                      /* Albedo */
         double *asymmetry                    /* Asymmetry */
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
