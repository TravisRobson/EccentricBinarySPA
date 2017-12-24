#ifndef Ecc_Adiabat_Evol_h
#define Ecc_Adiabat_Evol_h

#include "Ecc_IO.h"

int func(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);

void evolve_binary(struct EccBinary *eb, struct Data *data, double *y, gsl_odeiv2_driver *d, FILE *fptr);


int func_no_RR(double t, const double y[], double f[], void *params);
void evolve_no_RR(struct EccBinary *eb, struct Data *data, double y[], gsl_odeiv2_driver *d, FILE *fptr);
int jac_no_RR(double t, const double y[], double *dfdy, double dfdt[], void *params);
#endif /* Ecc_Adiabat_Evol_h */