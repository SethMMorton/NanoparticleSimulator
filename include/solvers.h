#ifndef SOLVERS_H
#define SOLVERS_H

#include <stdbool.h>
#include <complex.h>
#include "constants.h"

/* Main solver */
int npsolve (int nlayers,
             double rad[3],
             double rel_rad[nlayers][3],
             int indx[nlayers],
             double mrefrac,
             bool size_correct,
             bool cross_section,
             double extinct[NLAMBDA],
             double scat[NLAMBDA],
             double absorb[NLAMBDA]
           );

/* Quasistatic approx */
int quasi (int nlayers,
           double complex dielec[nlayers],
           double mdie,
           double rel_rad[nlayers][3],
           double rad[3],
           double size_param,
           double *extinct,
           double *scat,
           double *absorb
         );

/* Mie theory */
int mie (int nlayers,
         double complex refrac_indx[nlayers],
         double rel_rad[nlayers],
         double size_param,
         double *extinct,
         double *scat,
         double *absorb,
         double *backscat,
         double *rad_pressure,
         double *albedo,
         double *asymmetry
        );

#endif // SOLVERS_H
