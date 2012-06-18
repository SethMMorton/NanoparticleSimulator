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

#endif // SOLVERS_H
