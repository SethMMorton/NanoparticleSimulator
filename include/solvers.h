#ifndef SOLVERS_H
#define SOLVERS_H

/* Quasistatic approx */
int quasi (int nlayers,
           double complex dielec[nlayers],
           double mdie,
           double rel_rad[nlayers][3],
           double rad[3],
           double size_param,
           double *extinct,
           double *scat,
           double *absorb);

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
