#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "Ecc_IO.h"
#include "Detector.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage)
{
    double val = (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3.1f%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

void read_soln(double *time_series, FILE *fptr, struct Data *data)
{
	long i;
	long N, NFFT, left_pad;
	double t, phi, p, e, Forb;
		
	N        = data->N;
	NFFT     = data->NFFT;
	left_pad = data->left_pad;
	
	// initialize array to zeros
	for (i=0; i<2*NFFT; i++)
	{
		time_series[i] = 0.;
	}
	
	for (i=0; i<N; i++)
	{
		fscanf(fptr, "%lf%lf%lf%lf%lf%lf\n", &t, &phi, &p, &e, &Forb, &time_series[2*(i+left_pad)]);
	} 
		
	fclose(fptr);
	
	return;
}

void read_soln_no_RR(double *time_series, FILE *fptr, struct Data *data)
{
	long i;
	long N, NFFT, left_pad;
	double t, phi;
		
	N        = data->N;
	NFFT     = data->NFFT;
	left_pad = data->left_pad;
	
	// initialize array to zeros
	for (i=0; i<2*NFFT; i++)
	{
		time_series[i] = 0.;
	}
	
	for (i=0; i<N; i++)
	{
		fscanf(fptr, "%lf%lf%lf\n", &t, &phi, &time_series[2*(i+left_pad)]);
	} 
		
	fclose(fptr);
	
	return;
}

void print_dft(double *time_series, FILE *fptr, struct Data *data)
{
	long i;
	long NFFT, under_samp;
	long iRe, iIm;
	
	double f, Sn, dt;
	
	NFFT       = data->NFFT;
	under_samp = data->under_samp;
	dt         = data->dt;

	for (i=1; i<=NFFT/2; i++)
	{
		if (i%under_samp == 0)
		{
			iRe = 2*i;
			iIm = 2*i+1;
			f = (double)(i)/(double)NFFT/dt;
			
			Sn = get_Sn(f);
			
			// multiply by dt to give right units and size
			fprintf(fptr, "%e %e %e %e\n", f, time_series[iRe]*dt, time_series[iIm]*dt, Sn);
		}
	}
	fclose(fptr);
	
	return;
}

void print_spa(double *time_series, FILE *fptr, struct Data *data)
{
	long i;
	long NFFT, under_samp;
	long iRe, iIm;
	
	double f, Sn, dt;
	
	NFFT       = data->NFFT;
	under_samp = data->under_samp;
	dt         = data->dt;

	for (i=1; i<=NFFT/2; i++)
	{
		if (i%under_samp == 0)
		{
			iRe = 2*i;
			iIm = 2*i+1;
			f = (double)(i)/(double)NFFT/dt;
			
			Sn = get_Sn(f);
			
			fprintf(fptr, "%e %e %e %e\n", f, time_series[iRe], time_series[iIm], Sn);
		}
	}
	fclose(fptr);
	
	return;
}

void fill_num_series(double *num_series, struct Data *data)
{
	long i;
	
	FILE *in_file;
	
	in_file = fopen("soln.dat", "r");
	read_soln(num_series, in_file, data);
	
	gsl_fft_complex_radix2_forward(num_series, 1, data->NFFT);
	print_dft(num_series, fopen("dft.dat", "w"), data);
	
	for (i=0; i<2*data->NFFT; i++)
	{	// Correct the units to make a true FT
		num_series[i] *= data->dt;
	}
	
	return;
}