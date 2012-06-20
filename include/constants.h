#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Program parameters */
#define NLAMBDA 800
#define NMATERIAL 49

/* Physical Constants */
#define PI 3.14159265358979323846
#define HBAR 6.5821189916e-16

/* Conversions */
#define NM2EV(x) (1239.0 / x)
#define MPERS2EV(x, r) (x * HBAR / ( r * 1.0E-9 ))

/* Convenience */
#define SQR(x) ((x) * (x))
#define CMPLX(r, i) ((r) + I*(i))

#endif // CONSTANTS_H
