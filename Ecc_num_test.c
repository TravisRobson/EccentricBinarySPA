#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "Ecc_SPA.h"
#include "Ecc_Binary.h"
#include "Ecc_Adiabat_Evol.h"
#include "Constants.h"
#include "Detector.h"
#include "Ecc_IO.h"

void setup_EccBinary(struct EccBinary *eb);


int main(int argc, char *argv[])
{		
	long iRe, iIm;
	long i;

	double time_spent;
	double snr, Sn ,f;
	
	FILE *out_file;
	FILE *in_file;
	
	clock_t begin = clock();
	
	fprintf(stdout, "==============================================================\n\n");
	
	struct EccBinary *eb = malloc(sizeof(struct EccBinary));
	struct Data *data    = malloc(sizeof(struct Data));
	
	setup_EccBinary(eb);
	
	data->sr = 1024.;	  	        	 // number of samples per second
	data->dt = 1./data->sr; 	 	     // time between each samples
	data->Ny = 0.5/data->dt;		     // Nyquist frequency
	
	gsl_odeiv2_system sys = {func, jac, 3, eb};
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_bsimp, 1e-12, 1e-12, 0.);
	
	double y[3] = {0., eb->p0, eb->e0};
	
	out_file = fopen("soln.dat", "w");
	evolve_binary(eb, data, y, d, out_file);
	
	data->df = 1./(double)data->N/data->dt; // frequency resolution (true before zero padding)
	fprintf(stdout, "df: %f\n", data->df);
	
	fprintf(stdout, "\nN: %ld samples\n", data->N);
	
	gsl_odeiv2_driver_free (d);
	
	clock_t end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nEvolution and print runtime: %f sec\n", time_spent);
	begin = clock();
	
	data->NFFT     = (long)(pow(2, floor( log((double)(data->N))/log(2.)  ) + 1.));
	data->left_pad = (long)((data->NFFT - data->N)/2);
	fprintf(stdout, "FFT df: %lf\n", 1./data->NFFT/data->dt);
	fprintf(stdout, "TRU df: %lf\n", 1./data->N/data->dt);
	fprintf(stdout, "NFFT: %ld\n", data->NFFT);
	
	double *time_series; time_series = malloc(2*data->NFFT*sizeof(double));
	in_file = fopen("soln.dat", "r");
	read_soln(time_series, in_file, data);
	
	gsl_fft_complex_radix2_forward(time_series, 1, data->NFFT);
	
	data->under_samp = 1;
	
	out_file = fopen("ecc_dft.dat", "w");
	print_dft(time_series, out_file, data);

	snr = 0.;
	for (i=1; i<=data->NFFT/2; i++)
	{
		if (i%data->under_samp == 0)
		{
			iRe = 2*i;
			iIm = 2*i+1;
			f = (double)(i)/(double)data->NFFT/data->dt;
			
			Sn = get_Sn(f);
			
			snr += (time_series[iRe]*time_series[iRe] + time_series[iIm]*time_series[iIm])/Sn;
		}
	}
	snr *= 4.*data->dt/data->NFFT*data->under_samp;
	snr  = sqrt(snr); 
	
	fprintf(stdout, "SNR: %f\n", snr);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nDFT and print runtime: %f sec\n", time_spent);
	
	
	fprintf(stdout, "\n==============================================================\n");
	
	free(time_series);
	free(eb);
	free(data);

	return 0;
}

void setup_EccBinary(struct EccBinary *eb)
{
	// set angles associated with orbit orientation and sky location
	eb->beta  = 3.*M_PI/7.;
	eb->iota  = 3.*M_PI/7.;
	eb->theta = 3.*M_PI/7.;
	eb->phi   = 3.*M_PI/7.;
	eb->psi   = 3.*M_PI/7.;
	
	eb->c2beta = cos(2.*eb->beta);
	eb->s2beta = sin(2.*eb->beta);
	eb->ciota  = cos(eb->iota);
	eb->siota  = sin(eb->iota);
	
	// calculate the Antenna patterns for each polarization
	set_Antenna(eb);
	
	// set distance to source
	eb->R  = 410.*1.0e6*PC/C; // 410 Mpc
	eb->m1 = 10.*TSUN;
	eb->m2 = 10.*TSUN;
	set_m(eb);
	set_eta(eb);
	set_Mc(eb);
	set_mu(eb);
	
	eb->eLSO = 0.1e-1; // LSO eccentricity
	set_FLSO(eb);	   // LSO mean orbital frequency
	set_c0(eb);		   // constant of integration
	
	eb->F0 = 2; //0.0437;   //f_min/5.;
	fprintf(stdout, "F LSO: %.3f Hz\n", eb->FLSO);

	eb->e0 = solve_sigma_e(eb);
	fprintf(stdout, "e0: %.3f\n", eb->e0);
	
	// dimensionless semi-latus rectum
	eb->p0 = (1. - eb->e0*eb->e0)*pow(eb->m/PI2/PI2/eb->F0/eb->F0, 1./3.)/eb->m;
	
	// for the sake of the SPA
	eb->lc = 0.;
	eb->tc = 0.;

	return;
}




