#include "constants.h"

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


