#ifndef SOLVERS_H
#define SOLVERS_H

#include <stdbool.h>
#include <complex.h>
#include "constants.h"

/* Main solver */
int npsolve (int nlayers,
             double rad[XYZ],
             double rel_rad[MAXLAYERS][XYZ],
             int indx[MAXLAYERS],
             double mrefrac,
             bool size_correct,
             bool cross_section,
             double *sphere_rad,
             double extinct[NLAMBDA],
             double scat[NLAMBDA],
             double absorb[NLAMBDA]
           );

/* Quasistatic approx */
int quasi (int nlayers,
           double complex dielec[MAXLAYERS],
           double mdie,
           double rel_rad[MAXLAYERS][XYZ],
           double rad[XYZ],
           double size_param,
           double *extinct,
           double *scat,
           double *absorb
         );

/* Mie theory */
int mie (int nlayers,
         double complex refrac_indx[MAXLAYERS],
         double rel_rad[MAXLAYERS],
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
