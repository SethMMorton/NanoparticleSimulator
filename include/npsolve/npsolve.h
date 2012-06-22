#ifndef NPSOLVE_H
#define NPSOLVE_H

#include <stdbool.h>
#include "constants.h" /* This imports the constants  with npsolve.h */

/* Color matching data found in standard_color_matching.c */
extern const struct cie color_observers;
extern const double ICCMat[3][3];

/* Function to return the index of a material, from material_index.c */
int material_index(char *material);

/* Container to hold all spectra results */
typedef struct spectra_containers *Spectra;

/* The actual npsolve function declaration */
int npsolve (int nlayers,
             double rad[XYZ],
             double rel_rad[MAXLAYERS][XYZ],
             int indx[MAXLAYERS],
             double mrefrac,
             bool size_correct,
             Spectra spectra;
           );

#endif // NPSOLVE_H
