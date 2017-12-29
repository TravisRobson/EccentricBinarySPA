#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

#include "Ecc_SPA.h"
#include "Ecc_Binary.h"
#include "Ecc_Adiabat_Evol.h"
#include "Constants.h"
#include "Detector.h"
#include "Ecc_IO.h"
#include "Ecc_Math.h"

void setup_EccBinary(struct EccBinary *eb);
void construct_Data(struct EccBinary *eb, struct Data *data);
void fill_Fisher(struct EccBinary *eb, struct Data *data);


int main(int argc, char *argv[])
{		
	int l     = 0;
	int k_max = 100; // for lc search
	
	long i, j, k;
	
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
		eb->tc = -find_max_tc(num_series, spa_series, inv_ft, data);
		
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
		
		// increment counter for loading bar
		l++;
		percent = (double)l/(double)(k_max+1);
		printProgress(percent);

		// increment lc
		eb->lc += PI2/(double)k_max;	
	}

	eb->lc = lc_mm;
	eb->tc = tc_mm;
	spa_matrix_to_array(spa_mat, spa_series, eb, data);;
	//print_spa(spa_series, fopen("spa.dat", "w"), data);
	
	snr_spa = get_SNR(spa_series, data);
	fprintf(stdout, "\nFF: %f\n",   max_match);
	fprintf(stdout, "max tc: %e\n", tc_mm);
	fprintf(stdout, "max lc: %e\n", lc_mm);
	fprintf(stdout, "spa SNR: %f\n", snr_spa);

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nFF duration: %f sec\n", time_spent);
	
	eb->NP = 11;
	eb->params = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		eb->params[i] = 0.;
	}
	
	eb->params[0]  = log(eb->Mc/TSUN);
	eb->params[1]  = log(eb->F0/10.); 		  // 10. Hz	
	eb->params[2]  = log(eb->c0/1.);  	      // 1. Hz
	eb->params[3]  = eb->lc;
	eb->params[4]  = log(-eb->tc/10.); 		  // 10 sec
	eb->params[5]  = log(eb->R/(1.0e6*PC/C)); // 1 Mpc
	eb->params[6]  = eb->beta;
	eb->params[7]  = cos(eb->iota);
	eb->params[8]  = eb->phi;
	eb->params[9]  = cos(eb->theta);
	eb->params[10] = eb->psi;
	
	eb->Fisher = malloc(eb->NP*sizeof(double *));
	for (i=0; i<eb->NP; i++)
	{
		eb->Fisher[i] = malloc(eb->NP*sizeof(double));
	}
	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<eb->NP; k++)
		{
			eb->Fisher[i][k] = 0.;
		}	
	}
	
	struct EccBinary *eb_p = malloc(sizeof(struct EccBinary));
	struct EccBinary *eb_m = malloc(sizeof(struct EccBinary));
	setup_interp(eb_p);
	setup_interp(eb_m);
	
	eb_p->FLSO = eb->FLSO;
	eb_m->FLSO = eb->FLSO;
	
	eb_p->j_min = eb->j_min;
	eb_m->j_min = eb->j_min;
	eb_p->j_max = eb->j_max;
	eb_m->j_max = eb->j_max;
	
	eb_p->NP = 11;
	eb_m->NP = 11;
	eb_p->params = malloc(eb_p->NP*sizeof(double));
	eb_m->params = malloc(eb_m->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		eb_m->params[i] = 0.;
		eb_p->params[i] = 0.;
	}
	for (i=0; i<eb->NP; i++)
	{
		eb_p->params[i] = eb->params[i];
		eb_m->params[i] = eb->params[i];
	}
	map_array_to_params(eb_p);
	map_array_to_params(eb_m);
	set_Antenna(eb_p);
	set_Antenna(eb_m);
	
	double *spa_p = malloc(2*data->NFFT*sizeof(double));
	double *spa_m = malloc(2*data->NFFT*sizeof(double));
	
	double **spa_dif = malloc(eb->NP*sizeof(double *));
	for (i=0; i<eb->NP; i++) spa_dif[i] = malloc(2*data->NFFT*sizeof(double));
	
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_p[i]   = 0.;
		spa_m[i]   = 0.;
	}
	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<2*data->NFFT; k++)
		{
			spa_dif[i][k] = 0.;	
		}
	}
				
	double ep = 1.0e-6;
	long iRe, iIm; 

	data->under_samp = 10;
	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<eb->NP; k++)
		{
			eb_p->params[k] = eb->params[k];
			eb_m->params[k] = eb->params[k];
		}
		eb_p->params[i] += ep;
		eb_m->params[i] -= ep;
		map_array_to_params(eb_p);
		map_array_to_params(eb_m);
		set_Antenna(eb_p);
		set_Antenna(eb_m);
		
		fill_spa_series(spa_p, eb_p, data);
		fill_spa_series(spa_m, eb_m, data);
		//fprintf(stdout, "internal snr: %f\n", get_SNR(spa_m, data));
		
		for (k=1; k<=data->NFFT/2/data->under_samp; k++)
		{
			iRe = 2*(k*data->under_samp);
			iIm = 2*(k*data->under_samp)+1;
		
			spa_dif[i][iRe] = (spa_p[iRe] - spa_m[iRe])/(2.*ep);
			spa_dif[i][iIm] = (spa_p[iIm] - spa_m[iIm])/(2.*ep);
		}
	}
	free(spa_p);
	free(spa_m);
	free(eb_p->params);
	free(eb_m->params);
	free(eb_p);
	free(eb_m);
	
	for (i=0; i<eb->NP; i++)
	{
		for (k=i; k<eb->NP; k++)
		{
			eb->Fisher[i][k] = get_overlap(spa_dif[i], spa_dif[k], data);
		}
	}
	for (i=0; i<eb->NP; i++)
	{
		free(spa_dif[i]);
	}
	
	for (i=0; i<eb->NP; i++)
	{
		for (k=i; k<eb->NP; k++)
		{
			eb->Fisher[k][i] = eb->Fisher[i][k];
		}
	}	
	invert_matrix(eb->Fisher, eb->NP);
	
	fprintf(stdout, "\nFisher error estimates\n");
	fprintf(stdout, "------------------------------------\n");
	fprintf(stdout, "Mc   perc. error: %e\n", sqrt(eb->Fisher[0][0]));
	fprintf(stdout, "F0   perc. error: %e\n", sqrt(eb->Fisher[1][1]));
	fprintf(stdout, "c0   perc. error: %e\n", sqrt(eb->Fisher[2][2]));
	fprintf(stdout, "lc     abs error: %e\n", sqrt(eb->Fisher[3][3]));
	fprintf(stdout, "tc   perc. error: %e\n", sqrt(eb->Fisher[4][4]));
	fprintf(stdout, "R    perc. error: %e\n", sqrt(eb->Fisher[5][5]));
	fprintf(stdout, "beta   abs error: %e\n", sqrt(eb->Fisher[6][6]));
	fprintf(stdout, "ci     abs error: %e\n", sqrt(eb->Fisher[7][7]));
	fprintf(stdout, "phi    abs error: %e\n", sqrt(eb->Fisher[8][8]));
	fprintf(stdout, "ctheta abs error: %e\n", sqrt(eb->Fisher[9][9]));
	fprintf(stdout, "psi    abs error: %e\n", sqrt(eb->Fisher[10][10]));
	
	//fprintf(stdout, "\nvals: %e, %e, %e\n", eb->Fisher[5][5], eb->Fisher[8][8], eb->Fisher[9][9]);
	
	
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
	
	free(eb->params);
	for (i=0; i<eb->NP; i++)
	{
		free(eb->Fisher[i]);
	}

	return 0;
}

void fill_Fisher(struct EccBinary *eb, struct Data *data)
{

	return;
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
	eb->e0 = 0.8;
	eb->p0 = (1. - eb->e0*eb->e0)*pow(eb->m/PI2/PI2/eb->F0/eb->F0, 1./3.)/eb->m;
	eb->c0 = pow(calc_sigma(eb->e0), 3./2.)*eb->F0;
	eb->FLSO = pow((1. + eb->e0)/(6. + 2.*eb->e0), 3./2.)/PI2/eb->m;
 	
 	fprintf(stdout, "e0, p0, F0: %f %f %f\n", eb->e0, eb->p0, eb->F0);
	fprintf(stdout, "R: %e\n", eb->R);
	fprintf(stdout, "F LSO: %.3f Hz\n", eb->FLSO);

	eb->lc = 0.;
	eb->tc = 0.;
	
	eb->j_max = 100;
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






