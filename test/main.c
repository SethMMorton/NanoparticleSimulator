#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include "solvers.h"
#include "constants.h"

#define NLAYERS 2
#define XYZ 3
#define WAVELENGTH 500
#define MEDIUMDIE 1.0

// Declare Fortran routines
void quasif(int *n, double complex [*n], double*, double [(*n)*3], double [3], double*, double*, double*, double*);
void mief(int *n, double complex [*n], double [*n], double*, double*, double*, double*, double*, double*, double*, double*);

int main(int argc, char *argv[])
{
    
    // Define the system
    double rel_rad_spheroid[NLAYERS][XYZ] = { { 0.8, 0.8, 0.8 },
                                              { 0.2, 0.2, 0.2 }
                                            };
    double rel_rad_spheroid_f[XYZ*NLAYERS]= { 0.8, 0.2, 0.8, 0.2, 0.8, 0.2 };
    double rel_rad_sphere[NLAYERS] = { 0.8, 0.2 };
    double rad[XYZ] = { 10.0, 10.0, 10.0 };
    double size_param = 2.0 * PI * rad[0] / WAVELENGTH;

    // Define the dielectric constant, and the refractive index from that
    double complex dielec[NLAYERS] = { 1.2 + 0.5*I, 0.8 + 0.9*I };
    double complex refrac_indx[NLAYERS];
    for (int i = 0; i < NLAYERS; i++) {
        double tmp0 = cabs(dielec[i]);
        double tmp1 = sqrt(( tmp0 + creal(dielec[i]) ) / 2.0);
        double tmp2 = sqrt(( tmp0 - creal(dielec[i]) ) / 2.0);
        refrac_indx[i] = tmp1 + tmp2*I;
    }

    // Variables for passing to fortran
    int nlay = NLAYERS;
    double mdie = MEDIUMDIE;

    // The return variables
    double extinct, scat, abs, back, radpress, alb, asym;

    // Test the quasistatic approx
    int res = quasi(NLAYERS, dielec, MEDIUMDIE, rel_rad_spheroid, rad, size_param, &extinct, &scat, &abs);
    switch (res) {
        case 0:
            break;
        case 1:
            fprintf(stderr, "Too many layers for quasistatic approx\n");
            return 1;
            break;
        case 2:
            fprintf(stderr, "Axes 2 and 3 must be identical in quasi\n");
            return 1;
            break;
        default:
            fprintf(stderr, "Unknown error in quasi\n");
            return 1;
            break;
        
    }
    printf("Quasi C: Ext %.16f, Sca %.16f, Abs %.16f\n", extinct, scat, abs);
    extinct = 0.0; scat = 0.0; abs = 0.0;
    // Fortran version

    quasif(&nlay, dielec, &mdie, rel_rad_spheroid_f, rad, &size_param, &extinct, &scat, &abs);
    printf("Quasi F: Ext %.16f, Sca %.16f, Abs %.16f\n", extinct, scat, abs);
    extinct = 0.0; scat = 0.0; abs = 0.0;

    // Test Mie theory
    res = mie(NLAYERS, refrac_indx, rel_rad_sphere, size_param, &extinct, &scat, &abs, &back, &radpress, &alb, &asym);
    switch (res) {
        case 0:
            break;
        case 1:
            fprintf(stderr, "Product of size parameter and refractive index to large for at least one layer\n");
            return 1;
            break;
        default:
            fprintf(stderr, "Unknown error in mie\n");
            return 1;
            break;
    }
    printf("Mie C: Ext %.16f, Sca %.16f, Abs %.16f\n", extinct, scat, abs);
    extinct = 0.0; scat = 0.0; abs = 0.0;
    // Fortran version
    mief(&nlay, refrac_indx, rel_rad_sphere, &size_param, &extinct, &scat, &abs, &back, &radpress, &alb, &asym);
    printf("Mie F: Ext %.16f, Sca %.16f, Abs %.16f\n", extinct, scat, abs);

    // Main solver test
    int indx[NLAYERS] = { 0, 4 };
    res = npsolve(NLAYERS, rad, rel_rad_spheroid, indx, MEDIUMDIE, false, false, &extinct, &scat, &abs);
    printf("NPSolve: Ext %.16f, Sca %.16f, Abs %.16f\n", extinct, scat, abs);

    return 0;
}
