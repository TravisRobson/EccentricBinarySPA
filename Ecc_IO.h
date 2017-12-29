#ifndef Ecc_IO_h
#define Ecc_IO_h

struct Data
{
	long N, NFFT;
	long left_pad;
	long under_samp;
	
	double dt, df;
	double sr, Ny;

};
void read_soln_no_RR(double *time_series, FILE *fptr, struct Data *data);

void read_soln(double *time_series, FILE *fptr, struct Data *data);
void print_dft(double *time_series, FILE *fptr, struct Data *data);
void print_spa(double *time_series, FILE *fptr, struct Data *data);
void fill_num_series(double *num_series, struct Data *data);
void printProgress(double percentage);

#endif /* Ecc_IO */