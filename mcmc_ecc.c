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
	
	long i, k;
	
	double time_spent;
	double snr, snr_spa, match;

	double max_match = 0.;
	double tc_mm     = 0.;
	double lc_mm     = 0.;
	
	double *num_series;
	double *spa_series;
	
	clock_t begin = clock();
	clock_t end;
	
	fprintf(stdout, "==============================================================\n\n");
	
	struct EccBinary *eb = malloc(sizeof(struct EccBinary));
	struct Data *data    = malloc(sizeof(struct Data));
	
	setup_EccBinary(eb);
	setup_interp(eb);
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
	
	snr_spa = get_SNR(spa_series, data);
	fprintf(stdout, "\nFF: %f\n",   max_match);
	fprintf(stdout, "max tc: %e\n", tc_mm);
	fprintf(stdout, "max lc: %e\n", lc_mm);
	fprintf(stdout, "spa SNR: %f\n", snr_spa);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nFF duration: %f sec\n", time_spent);
	
	eb->NP = 11;//12;
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
	
	eb_p->NP = 11;//12;
	eb_m->NP = 11;// 12;
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

	data->under_samp = 64;
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

	
	////////////////////////	MCMC time	/////////////////////////////
	long m;
	long NMCMC = (long)1.0e3;

	double **history;
	history = malloc(NMCMC/10*sizeof(double *));
	for (i=0; i<NMCMC/10; i++)
	{
		history[i] = malloc(eb->NP*sizeof(double));
	}

	// CONSTRUCT EIGENBS FOR FISHER JUMPS
	double *evals, **evecs;
	evals = malloc(eb->NP*sizeof(double));
	evecs = malloc(eb->NP*sizeof(double *));
	for (i=0; i<eb->NP; i++) evecs[i] = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<eb->NP; k++) evecs[i][k] = 0.;
		evals[i] = 0.;
	}
    matrix_eigenstuff(eb->Fisher, evecs, evals, eb->NP); // fisher matrix gets inverted here I believe

	// SET UP SEEDS AND RANDOM NUMBER BS
	
	int seed = atoi(argv[1]);
	gsl_rng_env_setup();
	const gsl_rng_type *TT = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(TT);
	gsl_rng_set(r, seed);
	
	struct EccBinary *eb_x = malloc(sizeof(struct EccBinary));
	struct EccBinary *eb_y = malloc(sizeof(struct EccBinary));
	setup_interp(eb_x);
	setup_interp(eb_y);
	
	eb_x->FLSO = eb->FLSO;
	eb_y->FLSO = eb->FLSO;
	
	eb_x->j_min = eb->j_min;
	eb_y->j_min = eb->j_min;
	eb_x->j_max = eb->j_max;
	eb_y->j_max = eb->j_max;
	
	eb_x->NP = 11;//12;
	eb_y->NP = 11;//12;
	eb_x->params = malloc(eb->NP*sizeof(double));
	eb_y->params = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		eb_x->params[i] = 0.;
		eb_y->params[i] = 0.;
	}
	for (i=0; i<eb->NP; i++)
	{
		eb_x->params[i] = eb->params[i];
		eb_y->params[i] = eb->params[i];
	}
	map_array_to_params(eb_x);
	map_array_to_params(eb_y);
	set_Antenna(eb_x);
	set_Antenna(eb_y);
	
	FILE *out_file;
	double accept = 0.;
	out_file = fopen("chain.dat", "w");
	int meet_priors = 1;
	long itr = 0;
	
	struct EccBinary *eb_max = malloc(sizeof(struct EccBinary));
	setup_interp(eb_max);
	eb_max->FLSO = eb->FLSO;
	eb_max->j_min = eb->j_min;
	eb_max->j_max = eb->j_max;
	eb_max->NP = 11;//12;
	eb_max->params = malloc(eb->NP*sizeof(double));
	for (i=0; i<eb->NP; i++)
	{
		eb_max->params[i] = 0.;
	}
	for (i=0; i<eb->NP; i++)
	{
		eb_max->params[i] = eb->params[i];
	}
	map_array_to_params(eb_max);
	set_Antenna(eb_max);

	double logL_max = -10000.;
	double logLx, logLy, loga;
	
	double *spa_x, *spa_y;
	spa_x = malloc(2*data->NFFT*sizeof(double));
	spa_y = malloc(2*data->NFFT*sizeof(double));
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_x[i] = 0.;
		spa_y[i] = 0.;
	}
	
	double jump; int g, h; double alpha, beta;

 	logLx = get_logL(num_series, spa_series, data);
	fprintf(stdout, "\ninitial logLx: %f\n\n", logLx);
	
	for (m=0; m<NMCMC; m++)
	{
		if (m%(int)(NMCMC*0.01) == 0 && m!=0)
		{
			fprintf(stdout, "Percent Complete: %.0f %%		", (double)m/(double)NMCMC*100.);
			fprintf(stdout, "Acceptance Rate: %.1f %%\n",  accept/(double)m*100);
		} 
		
		jump = gsl_rng_uniform(r);
		if (jump < 0.5)
		{
			l = (int)(gsl_rng_uniform(r)*(double)eb->NP);
			for (k=0; k<eb->NP; k++)
			{
				eb_y->params[k] = eb_x->params[k] + gsl_ran_gaussian(r, 1.)*evecs[k][l]/sqrt(evals[l]*eb->NP);//*eb->NP/TT);
			}
		} else 
		{
			if (m/10<2) 
			{
				for (k=0; k<eb->NP; k++)
				{	// i.e. change nothing
					eb_y->params[k] = eb_x->params[k];
				}
			} else 
			{	
				g = (int)(((double)m/10*gsl_rng_uniform(r)));
				do {
					h = (int)(((double)m/10*gsl_rng_uniform(r)));
				} while (g==h);
				alpha = 1.0;
				beta = gsl_rng_uniform(r);
				if (beta < 0.9) alpha = gsl_ran_gaussian(r, 1.);
				for (k=0; k<eb->NP; k++)
				{
					eb_y->params[k] = eb_x->params[k] + alpha*(history[g][k] - history[h][k]);
				}
			}
		}
		if (eb_y->params[0]  <   log(0.1)    || eb_y->params[0]  > log(50.))   meet_priors = 0; // Mc
		if (eb_y->params[1]  <   log(0.001)  || eb_y->params[1]  > log(50.))   meet_priors = 0; // F0 
		if (eb_y->params[2]  <   log(0.001)  || eb_y->params[2]  > log(50.))   meet_priors = 0; // c0
		if (eb_y->params[3]  <=  0.          || eb_y->params[3]  > PI2)        meet_priors = 0; // lc
		if (eb_y->params[4]  <   log(1.0e-4) || eb_y->params[4]  > log(100.))  meet_priors = 0; // tc
		if (eb_y->params[5]  <   log(1.0e-6) || eb_y->params[5]  > log(1.0e6)) meet_priors = 0; // R
		if (eb_y->params[6]  <=  0.          || eb_y->params[6]  > PI2)        meet_priors = 0; // beta
		if (eb_y->params[7]  <= -1.          || eb_y->params[7]  > 1.)         meet_priors = 0; // ciota
		if (eb_y->params[8]  <=  0.          || eb_y->params[8]  > PI2)        meet_priors = 0; // phi
		if (eb_y->params[9]  <= -1.          || eb_y->params[9]  > 1.)         meet_priors = 0; // ctheta
		if (eb_y->params[10] <=  0.          || eb_y->params[10] > M_PI)       meet_priors = 0; // psi
	//	if (eb_y->params[11] <   log(0.001)  || eb_y->params[11] > log(50.))   meet_priors = 0; // FLSO

		//for (l=0; l<eb->NP; l++) fprintf(stdout, "y params[%d]: %f\n", l, eb_y->params[l]); 
		
		if (meet_priors == 1)
		{	// calculate signal associated with proposed source
			map_array_to_params(eb_y);
			set_Antenna(eb_y);
			fill_spa_series(spa_y, eb_y, data);
			logLy = get_logL(num_series, spa_y, data);	
			loga = log(gsl_rng_uniform(r));	
		}
		
		if (loga < (logLy - logLx) && meet_priors == 1)
		{
			accept += 1.;
			for (i=0; i<eb->NP; i++)
			{
				eb_x->params[i] = eb_y->params[i];
			}
			map_array_to_params(eb_x);
			set_Antenna(eb_x);
			logLx = logLy;
			if (logLx > logL_max)
			{
				logL_max = logLx;
				for (i=0; i<eb->NP; i++)
				{
					eb_max->params[i] = eb_x->params[i];
				}
				map_array_to_params(eb_max);
				set_Antenna(eb_max);
			}
		}

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
	
	fprintf(stdout, "\n[%ld]Acceptance rate: %.2f%%\n", m, 100.*accept/itr);
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
// 	free(spa_0);
	
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

	
	free(eb->params);
	for (i=0; i<eb->NP; i++)
	{
		free(eb->Fisher[i]);
	}
	
	for (i=0; i<NMCMC/10; i++)
	{
		free(history[i]);
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










