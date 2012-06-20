#ifndef MIE_HELPERS_H
#define MIE_HELPERS_H

#include <complex.h>

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

#endif // MIE_HELPERS_H
