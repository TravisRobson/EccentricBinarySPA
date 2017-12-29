#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "Ecc_SPA.h"
#include "Ecc_Binary.h"
#include "Ecc_Adiabat_Evol.h"
#include "Constants.h"
#include "Detector.h"
#include "Ecc_IO.h"
#include "Ecc_Math.h"

void setup_EccBinary(struct EccBinary *eb);
void construct_Data(struct EccBinary *eb, struct Data *data);


int main(int argc, char *argv[])
{		
	int l     = 0;
	int k_max = 100; // for lc search
	int m_max = 0;   // for tc search
	
	long i, j, k, m;
	
	double time_spent;
	double snr, snr_spa, match;

	double max_match = 0.;
	double tc_mm     = 0.;
	double lc_mm     = 0.;
	double percent   = 0.;
	
	double *num_series;
	double *inv_ft;
	double *spa_series;
	double *spa_0;
	double **spa_mat;
		
	clock_t begin = clock();
	clock_t end;
	
	fprintf(stdout, "==============================================================\n\n");
	
	struct EccBinary *eb = malloc(sizeof(struct EccBinary));
	struct Data *data    = malloc(sizeof(struct Data));
	
	setup_EccBinary(eb);
	setup_interp(eb);
	construct_Data(eb, data);
	
	num_series  = malloc(2*data->NFFT*sizeof(double));
	inv_ft      = malloc(2*data->NFFT*sizeof(double));
	spa_series  = malloc(2*data->NFFT*sizeof(double));
	spa_0		= malloc(2*data->NFFT*sizeof(double));
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_series[i] = 0.;
		spa_0[i]      = 0.;
	}
	//fill_spa_series(spa_0, eb, data);
	fill_num_series(num_series, data);
	snr = get_SNR(num_series, data);
	fprintf(stdout, "num SNR: %f\n\n", snr);
	
	begin = clock();
	spa_mat = malloc(eb->j_max*sizeof(double *));
	for (j=0; j<eb->j_max; j++)
	{
		spa_mat[j] = malloc(2*data->NFFT*sizeof(double));
	}
	for (j=0; j<eb->j_max; j++)
	{
		for (i=0; i<2*data->NFFT; i++)
		{
			spa_mat[j][i] = 0.;
		}
	}
	fill_SPA_matrix(spa_mat, eb, data);
	spa_matrix_to_array(spa_mat, spa_0, eb, data);
		
	for (k=0; k<=k_max; k++)
	{	
		eb->tc = 0.;	// reset tc	  
		spa_matrix_to_array(spa_mat, spa_series, eb, data); // Bc lc has been updated
		spa_matrix_to_array(spa_mat, spa_0, eb, data);		// Bc lc has been updated
		
		// iFFT to find tc which maximizes for current lc
		eb->tc  = -find_max_tc(num_series, spa_series, inv_ft, data);
		eb->tc -=  data->dt;
		
		// perform local search for tc max
		for (m=0; m<=m_max; m++)
		{
			fill_spa_series_new(spa_series, data, spa_0, eb);
			snr_spa = get_SNR(spa_series, data);
			match   = get_overlap(spa_series, num_series, data)/snr/snr_spa;
			
			// update maximum match
			if (fabs(match)>max_match)
			{
				max_match = fabs(match);
				tc_mm = eb->tc;
				lc_mm = eb->lc;
			} 
			
			// increment tc
			eb->tc += 2.*data->dt/(double)m_max;

			// increment counter for loading bar
			l++;
			percent = (double)l/(double)((k_max+1)*(m_max+1));
			printProgress(percent);
		}
		// increment lc
		eb->lc += PI2/(double)k_max;	
	}

	eb->lc = lc_mm;
	eb->tc = tc_mm;
	spa_matrix_to_array(spa_mat, spa_series, eb, data);;
	print_spa(spa_series, fopen("spa.dat", "w"), data);
	
	snr_spa = get_SNR(spa_series, data);
	fprintf(stdout, "\nFF: %f\n",    max_match);
	fprintf(stdout, "max tc: %e\n",  tc_mm);
	fprintf(stdout, "max lc: %e\n",  lc_mm);
	fprintf(stdout, "spa SNR: %f\n", snr_spa);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nFF duration: %f sec\n", time_spent);
	
	fprintf(stdout, "\n==============================================================\n");
	
	free(num_series);
	free(spa_series);
	free(spa_0);
	
	free(eb);
	free(data);
	for (i=0; i<eb->j_max; i++)
	{
		free(spa_mat[i]);
	}


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
	
	eb->eLSO = 0.01;   // LSO eccentricity

	eb->F0 = 3.0;
	eb->e0 = 0.7;
	eb->p0 = (1. - eb->e0*eb->e0)*pow(eb->m/PI2/PI2/eb->F0/eb->F0, 1./3.)/eb->m;
	eb->c0 = pow(calc_sigma(eb->e0), 3./2.)*eb->F0;
	eb->FLSO = pow((1. + eb->e0)/(6. + 2.*eb->e0), 3./2.)/PI2/eb->m;
 	
 	fprintf(stdout, "e0, p0: %e %e\n", eb->e0, eb->p0);
	fprintf(stdout, "R: %e\n", eb->R);
	fprintf(stdout, "F LSO: %.3f Hz\n", eb->FLSO);

	eb->lc = 0.;
	eb->tc = 0.;
	
	eb->j_max = 30;
	eb->j_min = 1;

	return;
}

void construct_Data(struct EccBinary *eb, struct Data *data)
{
	double time_spent;
	
	double y[3] = {0., eb->p0, eb->e0};
	
	clock_t begin = clock();
	
	FILE *out_file;
	
	data->sr = 8192;			  	     // number of samples per second
	data->dt = 1./data->sr; 	 	     // time between each samples
	data->Ny = 0.5*data->sr;		     // Nyquist frequency
	
	gsl_odeiv2_system sys = {func, jac, 3, eb};
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_bsimp, 1e-12, 1e-12, 0.);
	
	out_file = fopen("soln.dat", "w");
	evolve_binary(eb, data, y, d, out_file);
	
	// frequency resolution (true before zero padding)
	data->df = 1./(double)data->N/data->dt; 

	gsl_odeiv2_driver_free (d);
	
	clock_t end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	
	data->under_samp = 1;
	data->NFFT     = (long)(pow(2, floor( log((double)(data->N))/log(2.)  ) + 1.));
	data->left_pad = (long)((data->NFFT - data->N)/2);
	
	fprintf(stdout, "Evolution and print runtime: %f sec\n", time_spent);
	fprintf(stdout, "\nN: %ld samples\n", data->N);
	fprintf(stdout, "FFT df: %lf\n", 1./data->NFFT/data->dt);
	fprintf(stdout, "TRU df: %lf\n", 1./data->N/data->dt);
	fprintf(stdout, "NFFT: %ld\n", data->NFFT);

	return;
}





