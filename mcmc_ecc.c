#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Ecc_SPA.h"
#include "Ecc_Binary.h"
#include "Ecc_Adiabat_Evol.h"
#include "Ecc_IO.h"
#include "Ecc_Math.h"
#include "Constants.h"
#include "Detector.h"

void setup_EccBinary(struct EccBinary *eb);
void construct_Data(struct EccBinary *eb, struct Data *data);
void check_priors(struct EccBinary *eb_y, int *meet_priors);
void copy_EccBinary(struct EccBinary *dest, struct EccBinary *src);
void setup_xy_ebs(struct EccBinary *eb_x, struct EccBinary *eb_y, struct EccBinary *eb);
void fisher_jump(struct EccBinary *eb_x, struct EccBinary *eb_y, struct EccBinary *eb, gsl_rng *r);
void diff_ev_jump(struct EccBinary *eb_x, struct EccBinary *eb_y, struct EccBinary *eb, gsl_rng *r, double **history, long m);

int main(int argc, char *argv[])
{		
	int l = 0;
	
	long m;
	long i, k;
	long NMCMC = (long)1e4;
	
	double time_spent;
	double snr, snr_spa, match;
	double max_match;
	double jump; 
	double logL_max = -1.e30;
	double logLx, logLy, loga;
	
	double *num_series;
	double *spa_series;
	double **history;
	
	clock_t begin = clock();
	clock_t end;
	
	fprintf(stdout, "==============================================================\n\n");
	
	struct EccBinary *eb = malloc(sizeof(struct EccBinary));
	struct Data *data    = malloc(sizeof(struct Data));
	
	setup_EccBinary(eb);
	construct_Data(eb, data);
	
	num_series  = malloc(2*data->NFFT*sizeof(double));
	spa_series  = malloc(2*data->NFFT*sizeof(double));
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_series[i] = 0.;
	}
	fill_num_series(num_series, data);
	snr = get_SNR(num_series, data);
	fprintf(stdout, "num SNR: %f\n\n", snr);
	
	begin = clock();

	max_match = max_spa_tc_lc(eb, data, spa_series, num_series, snr);
	//print_spa(spa_series, fopen("spa.dat", "w"), data);
	eb->params[3]  = eb->lc;
	eb->params[4]  = log(-eb->tc/10.);
	snr_spa = get_SNR(spa_series, data);
	fprintf(stdout, "\nFF: %f\n",   max_match);
	fprintf(stdout, "max tc: %e\n", eb->tc);
	fprintf(stdout, "max lc: %e\n", eb->lc);
	fprintf(stdout, "spa SNR: %f\n", snr_spa);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nFF duration: %f sec\n", time_spent);
	
	data->under_samp = 64;
	calc_Fisher(eb, data);
		

	
	history = malloc(NMCMC/10*sizeof(double *));
	for (i=0; i<NMCMC/10; i++)
	{
		history[i] = malloc(eb->NP*sizeof(double));
	}

	eb->evals = malloc(eb->NP*sizeof(double));
	eb->evecs = malloc(eb->NP*sizeof(double *));
	for (i=0; i<eb->NP; i++) eb->evecs[i] = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<eb->NP; k++) eb->evecs[i][k] = 0.;
		eb->evals[i] = 0.;
	}
    matrix_eigenstuff(eb->Fisher, eb->evecs, eb->evals, eb->NP); // fisher matrix gets inverted here I believe

	int seed = atoi(argv[1]);
	gsl_rng_env_setup();
	const gsl_rng_type *TT = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(TT);
	gsl_rng_set(r, seed);
	
	struct EccBinary *eb_x = malloc(sizeof(struct EccBinary));
	struct EccBinary *eb_y = malloc(sizeof(struct EccBinary));
	setup_xy_ebs(eb_x, eb_y, eb);
	
	FILE *out_file;
	double accept = 0.;
	out_file = fopen("chain.dat", "w");
	int meet_priors = 1;
	long itr = 0;
	
	struct EccBinary *eb_max = malloc(sizeof(struct EccBinary));
	eb_max->NP = 11;//12;
	eb_max->params = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++) eb_max->params[i] = 0.;
	setup_interp(eb_max);
	copy_EccBinary(eb_max, eb);

	
	double *spa_x, *spa_y;
	spa_x = malloc(2*data->NFFT*sizeof(double));
	spa_y = malloc(2*data->NFFT*sizeof(double));
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_x[i] = 0.;
		spa_y[i] = 0.;
	}

 	logLx = get_logL(num_series, spa_series, data);
 	logLy = -1.0e30;
 	loga  = 0.;
	fprintf(stdout, "\ninitial logLx: %f\n\n", logLx);
	fprintf(stdout, "MCMC Initialized\n\n");
	for (m=0; m<NMCMC; m++)
	{
		if (m%(int)(NMCMC*0.01) == 0 && m!=0)
		{
			printProgress((double)m/(double)(NMCMC));
		} 
		
		jump = gsl_rng_uniform(r);
		if (jump < 0.5) 
		{
			fisher_jump(eb_x, eb_y, eb, r);
		} else 
		{
			diff_ev_jump(eb_x, eb_y, eb, r, history, m);
		}

		check_priors(eb_y, &meet_priors);
		//fprintf(stdout, "meet_priors[%ld]: %d\n", m/10, meet_priors);
		if (meet_priors == 1)
		{	// calculate signal associated with proposed source
			map_array_to_params(eb_y);
			set_Antenna(eb_y);
			fill_spa_series(spa_y, eb_y, data);
			
			logLy = get_logL(num_series, spa_y, data);	
			loga  = log(gsl_rng_uniform(r));	
		}
		//fprintf(stdout, "logLy[%ld]: %f\n", m, logLy);
		if (logLy > -INFINITY){//fprintf(stdout, "logLy[%ld]: %f\n", m, logLy);
		if (loga < (logLy - logLx) && meet_priors == 1)
		{
			accept += 1.;
			copy_EccBinary(eb_x, eb_y);
			logLx = logLy;
			
			if (logLx > logL_max)
			{
				logL_max = logLx;
				copy_EccBinary(eb_max, eb_x);
			}
		}}

		if (itr%10 == 0)
		{
			fprintf(out_file, "%ld %.12g ", itr/10, logLx);
			for (l=0; l<eb->NP; l++)
			{
				fprintf(out_file, "%.12g ", eb_x->params[l]);
				history[itr/10][l] = eb_x->params[l];
			}
			fprintf(out_file, "\n");
			
		}
		itr += 1;
		meet_priors = 1; // reset
	}
	
	fclose(out_file);
	printProgress((double)m/(double)(NMCMC));
	
	fprintf(stdout, "\n\nAcceptance rate: %.2f%%\n", 100.*accept/itr);
	fprintf(stdout, "logL_max: %f\n", logL_max);
	
	double *spa_max;
	spa_max = malloc(2*data->NFFT*sizeof(double));
	fill_spa_series(spa_max, eb_max, data);
	snr_spa = get_SNR(spa_max, data); snr = get_SNR(num_series, data); // repeat bc change in under sampling
	match   = get_overlap(spa_max, num_series, data)/snr/snr_spa;
	
	fprintf(stdout, "FF: %f\n", match);
	
	fprintf(stdout, "Mc   dif: %e\n", exp(eb->params[0])*TSUN  - exp(eb_max->params[0])*TSUN);
	fprintf(stdout, "F0   dif: %e\n", exp(eb->params[1])*10.   - exp(eb_max->params[1])*10.);
	fprintf(stdout, "c0   dif: %e\n", exp(eb->params[2])*1.    - exp(eb_max->params[2])*1.);
//	fprintf(stdout, "FLSO dif: %e\n", exp(eb->params[11])*100. - exp(eb_max->params[11])*100.);
		
	
	fprintf(stdout, "\n==============================================================\n");
	
	free(spa_max);
	free(num_series);
	free(spa_series);
	
	free(eb->params);
	for (i=0; i<eb->NP; i++) free(eb->Fisher[i]);
	free(eb);
	free(data);
	
	free(eb_x->params);
	free(eb_y->params);
	free(eb_max->params);
	free(eb_x);
	free(eb_y);
	free(eb_max);
	
	free(spa_x);
	free(spa_y);
	
	for (i=0; i<NMCMC/10; i++) free(history[i]);

	return 0;
}

void setup_EccBinary(struct EccBinary *eb)
{
	long i;
	
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
 	
 	fprintf(stdout, "e0, p0, F0: %f %f %f\n", eb->e0, eb->p0, eb->F0);
	fprintf(stdout, "R: %e\n", eb->R);
	fprintf(stdout, "F LSO: %.3f Hz\n", eb->FLSO);

	eb->lc = 0.;
	eb->tc = 0.;
	
	eb->j_max = 15;
	eb->j_min = 1;
	
	// set up parameter array
	eb->NP = 11; //12;
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
	//eb->params[11] = log(eb->FLSO/100.);
	
	setup_interp(eb);
	
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

void check_priors(struct EccBinary *eb_y, int *meet_priors)
{
	if (eb_y->params[0]  <   log(0.1)    || eb_y->params[0]  > log(50.))   *meet_priors = 0; // Mc
	if (eb_y->params[1]  <   log(0.001)  || eb_y->params[1]  > log(50.))   *meet_priors = 0; // F0 
	if (eb_y->params[2]  <   log(0.001)  || eb_y->params[2]  > log(50.))   *meet_priors = 0; // c0
	if (eb_y->params[3]  <=  0.          || eb_y->params[3]  > PI2)        *meet_priors = 0; // lc
	if (eb_y->params[4]  <   log(1.0e-4) || eb_y->params[4]  > log(100.))  *meet_priors = 0; // tc
	if (eb_y->params[5]  <   log(1.0e-6) || eb_y->params[5]  > log(1.0e6)) *meet_priors = 0; // R
	if (eb_y->params[6]  <=  0.          || eb_y->params[6]  > PI2)        *meet_priors = 0; // beta
	if (eb_y->params[7]  <= -1.          || eb_y->params[7]  > 1.)         *meet_priors = 0; // ciota
	if (eb_y->params[8]  <=  0.          || eb_y->params[8]  > PI2)        *meet_priors = 0; // phi
	if (eb_y->params[9]  <= -1.          || eb_y->params[9]  > 1.)         *meet_priors = 0; // ctheta
	if (eb_y->params[10] <=  0.          || eb_y->params[10] > M_PI)       *meet_priors = 0; // psi
	
//	if (eb_y->params[11] <   log(0.001)  || eb_y->params[11] > log(50.))   meet_priors = 0; // FLSO
	
	return;
}

void copy_EccBinary(struct EccBinary *dest, struct EccBinary *src)
{
	int i;
	 
	dest->FLSO  = src->FLSO;
	dest->j_min = src->j_min;
	dest->j_max = src->j_max;
	dest->NP    = src->NP;

	for (i=0; i<src->NP; i++)
	{
		dest->params[i] = src->params[i];
	}
	map_array_to_params(dest);
	set_Antenna(dest);
	
	return;
}

void setup_xy_ebs(struct EccBinary *eb_x, struct EccBinary *eb_y, struct EccBinary *eb)
{	
	int i;
	
	setup_interp(eb_x);
	setup_interp(eb_y);
	eb_x->params = malloc(eb->NP*sizeof(double));
	eb_y->params = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		eb_x->params[i] = 0.;
		eb_y->params[i] = 0.;
	}
	
	copy_EccBinary(eb_x, eb);
	copy_EccBinary(eb_y, eb);
	
	return;
}

void fisher_jump(struct EccBinary *eb_x, struct EccBinary *eb_y, struct EccBinary *eb, gsl_rng *r)
{	
	int i, j;
	
	double *jump;
	jump = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++) jump[i] = 0.;
	
	j = (int)(gsl_rng_uniform(r)*(double)eb_y->NP);
	for (i=0; i<eb->NP; i++)
	{
		jump[i] = gsl_ran_gaussian(r, 1.)*eb->evecs[i][j]/sqrt(eb->evals[j]*eb->NP);
	}
	
	
	for (i=0; i<eb->NP; i++)
	{
		eb_y->params[i] = eb_x->params[i] + jump[i];
	}
	
	return;
}

void diff_ev_jump(struct EccBinary *eb_x, struct EccBinary *eb_y, struct EccBinary *eb, gsl_rng *r, double **history, long m)
{
	int i, j, k;
		
	double alpha, beta;
	
	if (m/10<2) 
	{
		for (i=0; i<eb->NP; i++)
		{	// i.e. change nothing
			eb_y->params[i] = eb_x->params[i];
		}
	} else 
	{	
		j = (int)(((double)m/10*gsl_rng_uniform(r)));
		do {
			k = (int)(((double)m/10*gsl_rng_uniform(r)));
		} while (j==k);
		
		alpha = 1.0;
		beta = gsl_rng_uniform(r);
		if (beta < 0.9) alpha = gsl_ran_gaussian(r, 1.);
		
		for (i=0; i<eb->NP; i++)
		{
			eb_y->params[i] = eb_x->params[i] + alpha*(history[j][i] - history[k][i]);
		}
	}
	
	return;
}


