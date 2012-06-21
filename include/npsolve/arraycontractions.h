#ifndef ARRAYCONTRACTIONS_H
#define ARRAYCONTRACTIONS_H

#include <complex.h>

/* Summation functions */
double sum(int n, double array[n]);
double sumrange(int n, double array[n], int start, int end);
double sum2(int m, int n, double array[m][n]);
double sum2range(int m, int n, double array[m][n], 
                 int mstart, int mend, int nstart, int nend);
double complex csum(int n, double complex array[n]);
double complex csumrange(int n, double complex array[n], int start, int end);
double complex csum2(int m, int n, double complex array[m][n]);
double complex csum2range(int m, int n, double complex array[m][n], 
                          int mstart, int mend, int nstart, int nend);

/* Product functions */
double prod(int n, double array[n]);
double prodrange(int n, double array[n], int start, int end);
double prod2(int m, int n, double array[m][n]);
double prod2range(int m, int n, double array[m][n], 
                  int mstart, int mend, int nstart, int nend);
double complex cprod(int n, double complex array[n]);
double complex cprodrange(int n, double complex array[n], int start, int end);
double complex cprod2(int m, int n, double complex array[m][n]);
double complex cprod2range(int m, int n, double complex array[m][n], 
                           int mstart, int mend, int nstart, int nend);


#endif // ARRAYCONTRACTIONS_H
