#ifndef Ecc_SPA_h
#define Ecc_SPA_h

#include "Ecc_Binary.h"


// calculate the mean anomaly Fourier series coefficients
double get_Cpj(struct EccBinary *eb, double j, double e);
double get_Spj(struct EccBinary *eb, double j, double e);
double get_Ccj(struct EccBinary *eb, double j, double e);
double get_Scj(struct EccBinary *eb, double j, double e);

// invert f(e;j;c0) relation, i.e. get stationary eccentricity
double get_es(struct EccBinary *eb, int j, double f);

// calc It(e) and Il(e) (with Hypergeometric functions or our approximation)
double get_It(struct EccBinary *eb, double e);

double get_Il(struct EccBinary *eb, double e);

// get phase for a given e, f and j
double get_phase(struct EccBinary *eb, int j, double e, double f);

// calculate the amplitude factor which is j dependent
double get_Aj(struct EccBinary *eb, double j, double e);

// get a harmonic component
void get_harmonic(struct EccBinary *eb, int j, double f, double *hjRe, double *hjIm);

// calculate the real an imaginary parts of the SPA at a given f
void get_eccSPA(struct EccBinary *eb, double f, double *spaRe, double *spaIm);


double term(int r, double x);

void setup_interp(struct EccBinary *eb);
double interp_phase(struct EccBinary *eb, double j, double alpha, double f);
double interp_es(struct EccBinary *eb, double alpha);

#endif /* Ecc_SPA_h */