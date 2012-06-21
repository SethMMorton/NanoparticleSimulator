/*=========== Light Scattering by Spherical Particles ==============
*   Calculations of extinction, scattering, absorption, etc.
*   efficiency factors for n-layered spheres.
*------------------------------------------------------------------
*  Input data  :
*   nlayers    : number of layers
*   rad        : the radius of each layer (from inner surface
*                                          to outer surface)
*   rel_rad    : The relative 
*   lambda     : the wavelength limits to scan over
*   nlambda    : number of wavelengths to calculate between limits
*   mrefrac    : refractive index of surrounding medium (default = 1)
*   lexact     : use exact theory to calculate (default = true)
*   drude      : the drude function parameters (optional)
*                1 = bound dielectric
*                2 = plasma frequency
*                3 = fermi velocity
*                4 = inverse lifetime (gamma)
*   nexp       : number of experimental wavelengths (optional)
*   expdie     : the experimental dielectric constant (optional)
*                At least one of drude and expdie must be given.
*                If both expdie and drude are given, drude is used
*                as a correction term for expdie.
*                1 = wavelength (nm)
*                2 = real dielelectric
*                3 = imaginary dielelectric
*                4 = real second derivative (for cubic interpolation)
*                5 = imaginary second derivative
*------------------------------------------------------------------
*  Output data :
*     wavelengths  : the wavelengths calculated
*     extinct  : extinction efficiency (Qext)
*     scat     : scattering efficiency (Qsca)
*     absorb   : absorption efficiency (Qabs)
*------------------------------------------------------------------
* Recursive algorithms of Wu & Wang (Radio Sci. 26, 1393, 1991)
* created by N.V. Voshchinnikov
* (c) 1999 Sobolev Astronomical Institute, St. Petersburg Univ.
*
* Adapted for C by Seth Morton, 2012.
* Added quasistatic approximation
*==================================================================*/

#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include "npsolve/constants.h"
#include "npsolve/arraycontractions.h"
#include "npsolve/wavelengths.h"
#include "npsolve/experimental_dielectrics.h"
#include "npsolve/drude_parameters.h"
#include "npsolve/solvers.h"
#include "npsolve/npsolve.h"

/* Container to hold the three spectra types */
typedef struct {
    double extinct[NLAMBDA];
    double scat[NLAMBDA];
    double absorb[NLAMBDA];
} SpectraTypes;

/* Container to hold all spectra results */
struct spectra_containers {
    SpectraTypes efficiencies;
    SpectraTypes cross_sections;
    SpectraTypes absorptivity;
    SpectraTypes solution;
};

/* Constant variables from the headers */
extern const double wavelengths[NLAMBDA];
extern const double drude_parameters[NMATERIAL][3];
extern const double complex experimental_dielectrics[NMATERIAL][NLAMBDA];

/* Drude function macro */
#define DRUDE(x, p, g, s) (p*p / ( x * ( x + I * ( g + s ) ) ))

int npsolve (int nlayers,                    /* Number of layers */
             double rad[XYZ],                /* Radius of object */
             double rel_rad[MAXLAYERS][XYZ], /* Relative radii of layers */
             int indx[MAXLAYERS],            /* Material index of layers */
             double mrefrac,                 /* Refractive index of medium */
             bool size_correct,              /* Use size correction? */
             Spectra spectra;                /* Results */
           )
{

    /* Dielectric function and refractive index */
    double complex dielec[nlayers], refrac_indx[nlayers];

    /* If the second or third components are negative, it is a sphere
       and thus Mie theory is used.  Otherwise, quasistatic is used. */
    bool lmie = false;
    if (rad[1] < 0.0 || rad[2] < 0.0)
        lmie = true;

    /* Calculate radius of a sphere with an equivalent volume */
    double *sphere_rad;
    if (lmie)
        *sphere_rad = rad[0];
    else
        *sphere_rad = pow(prod(3, rad), ( 1.0 / 3.0 ));

    /***************************************************
     * Loop over each wavelength to calculate properties
     ***************************************************/

    for (int i = 0; i < NLAMBDA; i++) {

        /* Determine size parameter */
        double size_param = 2.0 * PI * (*sphere_rad) * mrefrac / wavelengths[i];

        /* Skip if size_param is too small */
        if (size_param < 0.1E-6)
            continue;

        /*****************************************************************
         * Calculate dielectric constant & refractive index for each layer
         *****************************************************************/

        for (int j = 0; j < nlayers; j++) {

            /* Grab dielectric from experiment */
            dielec[j] = experimental_dielectrics[indx[j]][i];

            /* Correct for size if the asked for */
            if (size_correct) {

                /* Extract the drude parameters */
                double pf = drude_parameters[indx[j]][0];
                double gm = drude_parameters[indx[j]][1];
                double sc = drude_parameters[indx[j]][2];

                /* Energy in electron volts (omega) */
                double om = NM2EV(wavelengths[i]);

                /* Use the drude model to size-correct experimental data */
                dielec[j] = dielec[j]
                          - DRUDE(om, pf, gm, 0.0)
                          + DRUDE(om, pf, gm, MPERS2EV(sc, *sphere_rad));

            }

            /* Turn dielectric into refractive index if Mie theory */
            if (lmie) {
                double tmp0 = cabs(dielec[j]);
                double tmp1 = sqrt(( tmp0 + creal(dielec[j]) ) / 2.0);
                double tmp2 = sqrt(( tmp0 - creal(dielec[j]) ) / 2.0);
                refrac_indx[j] = CMPLX(tmp1, tmp2);
            }

        }

        /* Solve using the appropriate inputs and theory */
        if (lmie) {
            /* Relative radius for sphere */
            double srrad[nlayers];
            for (int i = 0; i < nlayers; i++)
                srrad[i] = rel_rad[i][0];
            double backscat, rad_pressure, albedo, asymmetry;
            int retval = mie(nlayers, refrac_indx, srrad, size_param,
                             &extinct[i], &scat[i], &absorb[i],
                             &backscat, &rad_pressure, &albedo, &asymmetry);
            if (retval == 1)
                return 1; /* Product of size param & ref index too large */
        } else {
            int retval = quasi(nlayers, dielec, SQR(mrefrac), rel_rad,
                               rad, size_param,
                               &extinct[i], &scat[i], &absorb[i]);
            if (retval == 1)
                return 2; /* Too many layers for quasistatic approx */
            else if (retval == 2)
                return 3; /* Axes 2 and 3 not identical */
        }
        /* Calculate cross sections */
        if (cross_section) {
            extinct[i] *= PI * SQR(*sphere_rad);
            scat[i]    *= PI * SQR(*sphere_rad);
            absorb[i]  *= PI * SQR(*sphere_rad);
        }

    }

    return 0;

}
