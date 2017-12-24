#ifndef Ecc_Math_h
#define Ecc_Math_h

#include "Ecc_Binary.h"
#include "Ecc_IO.h"

double get_SNR(double *ft, struct Data *data);
double get_overlap(double *a, double *b, struct Data *data);

#endif /* Ecc_Math_h */