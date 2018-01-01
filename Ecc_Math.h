#ifndef Ecc_Math_h
#define Ecc_Math_h

#include "Ecc_Binary.h"
#include "Ecc_IO.h"

double get_SNR(double *ft, struct Data *data);
double get_overlap(double *a, double *b, struct Data *data);
double get_logL(double *a, double *b, struct Data *data);
double find_max_tc(double *a, double *b, double *inv_ft, struct Data *data);
double max_spa_tc_lc(struct EccBinary *eb, struct Data *data, double *spa_series, double *num_series, double snr);

void fill_spa_series_new(double *spa_series, struct Data *data, double *spa_0, struct EccBinary *eb);
void fill_SPA_matrix(double **spa_matrix, struct EccBinary *eb, struct Data *data);
void spa_matrix_to_array(double **spa_matrix, double *spa_series, struct EccBinary *eb, struct Data *data);
void fill_spa_series(double *spa_series, struct EccBinary *eb, struct Data *data);

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N);
double invert_matrix(double **matrix, int N);


#endif /* Ecc_Math_h */