#include <stdlib.h>
#include <stdio.h>
#include "arraycontractions.h"

/***************************
 * Functions to sum an array
 ***************************/

/* Sum doubles */
double sum(int n, double array[n]) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s += array[i];
    }
    return s;
}

/* Sum doubles in range */
double sumrange(int n, double array[n], int start, int end) {
    double s = 0.0;
    for (int i = start; i < end + 1; i++) {
        s += array[i];
    }
    return s;
}

/* Sum doubles in 2-D */
double sum2(int m, int n, double array[m][n]) {
    double s = 0.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            s += array[i][j];
        }
    }
    return s;
}

/* Sum doubles in 2-D in range*/
double sum2range(int m, int n, double array[m][n], 
                 int mstart, int mend, int nstart, int nend) {
    double s = 0.0;
    for (int i = mstart; i < mend + 1; i++) {
        for (int j = nstart; j < nend + 1; j++) {
            s += array[i][j];
        }
    }
    return s;
}

/* Sum double complexs */
double complex csum(int n, double complex array[n]) {
    double complex s = 0.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        s += array[i];
    }
    return s;
}

/* Sum double complexs in range */
double complex csumrange(int n, double complex array[n], int start, int end) {
    double complex s = 0.0 + 0.0 * I;
    for (int i = start; i < end + 1; i++) {
        s += array[i];
    }
    return s;
}

/* Sum double complexs in 2-D */
double complex csum2(int m, int n, double complex array[m][n]) {
    double complex s = 0.0 + 0.0 * I;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            s += array[i][j];
        }
    }
    return s;
}

/* Sum double complexs in 2-D in range*/
double complex csum2range(int m, int n, double complex array[m][n], 
                          int mstart, int mend, int nstart, int nend) {
    double complex s = 0.0 + 0.0 * I;
    for (int i = mstart; i < mend + 1; i++) {
        for (int j = nstart; j < nend + 1; j++) {
            s += array[i][j];
        }
    }
    return s;
}

/*******************************************
 * Functions to take the product of an array
 *******************************************/

/* Multiply doubles */
double prod(int n, double array[n]) {
    double s = 1.0;
    for (int i = 0; i < n; i++) {
        s *= array[i];
    }
    return s;
}

/* Multiply doubles in range */
double prodrange(int n, double array[n], int start, int end) {
    double s = 1.0;
    for (int i = start; i < end + 1; i++) {
        s *= array[i];
    }
    return s;
}

/* Multiply doubles in 2-D */
double prod2(int m, int n, double array[m][n]) {
    double s = 1.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            s *= array[i][j];
        }
    }
    return s;
}

/* Multiply doubles in 2-D in range*/
double prod2range(int m, int n, double array[m][n], 
                  int mstart, int mend, int nstart, int nend) {
    double s = 1.0;
    for (int i = mstart; i < mend + 1; i++) {
        for (int j = nstart; j < nend + 1; j++) {
            s *= array[i][j];
        }
    }
    return s;
}

/* Multiply double complexs */
double complex cprod(int n, double complex array[n]) {
    double complex s = 1.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        s *= array[i];
    }
    return s;
}

/* Multiply double complexs in range */
double complex cprodrange(int n, double complex array[n], int start, int end) {
    double complex s = 1.0 + 0.0 * I;
    for (int i = start; i < end + 1; i++) {
        s *= array[i];
    }
    return s;
}

/* Multiply double complexs in 2-D */
double complex cprod2(int m, int n, double complex array[m][n]) {
    double complex s = 1.0 + 0.0 * I;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            s *= array[i][j];
        }
    }
    return s;
}

/* Multiply double complexs in 2-D in range*/
double complex cprod2range(int m, int n, double complex array[m][n], 
                    int mstart, int mend, int nstart, int nend) {
    double complex s = 1.0 + 0.0 * I;
    for (int i = mstart; i < mend + 1; i++) {
        for (int j = nstart; j < nend + 1; j++) {
            s *= array[i][j];
        }
    }
    return s;
}
