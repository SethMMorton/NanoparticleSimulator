#include <stdio.h>
#include <complex.h>
#include "solvers.h"

int main(int argc, char *argv[])
{
    
    double extinct, scat, abs;
    double complex dielec[2] = {1.2 + 0.5*I, 0.8 + 0.9*I};
    double rel_rad[2][3];
    rel_rad[0][0] = 0.8;
    rel_rad[0][1] = 0.8;
    rel_rad[0][2] = 0.8;
    rel_rad[1][0] = 1.0;
    rel_rad[1][1] = 1.0;
    rel_rad[1][2] = 1.0;
    double rad[3];
    rad[0] = 15.0;
    rad[1] = 10.0;
    rad[2] = 10.0;
    int res = quasi(2, dielec, 1.0, rel_rad, rad, 1.5, &extinct, &scat, &abs);
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
    printf("Ext %.16f, Sca %.16f, Abs %.16f\n", extinct, scat, abs);

    return 0;
}
