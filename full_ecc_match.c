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


int main(int argc, char *argv[])
{		
	long iRe, iIm;
	long i, k, m;

	double time_spent;
	double snr, Sn ,f;
	double snr_spa;
	double max_corr;
	double tc_max;
	
	double *time_series;
	double *inv_ft;
	double *spa_series; 
	
	FILE *in_file;
	
	clock_t begin = clock();
	clock_t end;
	
	fprintf(stdout, "==============================================================\n\n");
	
	struct EccBinary *eb = malloc(sizeof(struct EccBinary));
	struct Data *data    = malloc(sizeof(struct Data));
	
	
	
	setup_EccBinary(eb);
	setup_interp(eb);
	construct_Data(eb, data);
	
	time_series = malloc(2*data->NFFT*sizeof(double));
	inv_ft      = malloc(2*data->NFFT*sizeof(double));
	spa_series = malloc(2*data->NFFT*sizeof(double));
	
	in_file = fopen("soln.dat", "r");
	read_soln(time_series, in_file, data);
	gsl_fft_complex_radix2_forward(time_series, 1, data->NFFT);
	print_dft(time_series, fopen("dft.dat", "w"), data);
	for (i=0; i<2*data->NFFT; i++)
	{	// Correct the units to make a true FT
		time_series[i] *= data->dt;
	}
	
	snr = get_SNR(time_series, data);
	fprintf(stdout, "num SNR: %f\n\n", snr);

	
	begin = clock();
	double *harm_t_ser; harm_t_ser = malloc(2.*data->N*sizeof(double));
	for (i=0; i<data->N;i++)
	{
		harm_t_ser[i] = 0.;
	}
	double match = 0.;

	 
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_series[i] = 0.;
	}

	long j_max, j_min;
	j_max = 30; j_min = 1;
	double spaRe, spaIm;
	for (i=1; i<=data->NFFT/2; i++)
	{
		if (i%data->under_samp == 0)
		{
			f = (double)(i)/(double)data->NFFT/data->dt;
			get_eccSPA(eb, f, &spaRe, &spaIm, j_max, j_min);
			
			iRe = 2*i;
			iIm = 2*i+1;

			spa_series[iRe] = spaRe;
			spa_series[iIm] = spaIm;
		}
	}
	snr_spa = get_SNR(spa_series, data);
	
	for (i=1; i<=data->NFFT/2; i++)
	{
		if (i%data->under_samp == 0)
		{
			iRe = 2*i;
			iIm = 2*i+1;
			f = (double)(i)/(double)data->NFFT/data->dt;
		
			Sn = get_Sn(f);
		
			match += (spa_series[iRe]*time_series[iRe] + spa_series[iIm]*time_series[iIm])/Sn;
		}
	}
	match *= 4./(double)data->NFFT*(double)data->under_samp/data->dt;
	match /= snr*snr_spa;

	long m_max = 50;
	double max_match = 0.; double tc_mm = 0.; double lc_mm = 0.;
	int k_max = 50; long l = 0; double percent = 0.;
	
	for (k=1; k<=k_max+1; k++)
	{	
		// sample the SPA at current tc and lc
		eb->tc = 0.;
		for (i=1; i<=data->NFFT/2; i++)
		{
			if (i%data->under_samp == 0)
			{
				f = (double)(i)/(double)data->NFFT/data->dt;
				get_eccSPA(eb, f, &spaRe, &spaIm, j_max, j_min);
			
				iRe = 2*i;
				iIm = 2*i+1;

				spa_series[iRe] = spaRe;
				spa_series[iIm] = spaIm;
			}
		}
		
		// iFFT to find tc which maximizes for current lc
		max_corr = 0.;
		tc_max   = 0.;
		for (i=0; i<2*data->NFFT; i++)
		{
			inv_ft[i] = 0.;
		}
		for (i=1; i<=data->NFFT/2; i++)
		{
			f = (double)(i)/(double)data->NFFT/data->dt;
			Sn = get_Sn(f);
		
			iRe = 2*i;
			iIm = 2*i+1;

			inv_ft[iRe] = (time_series[iRe]*spa_series[iRe] + time_series[iIm]*spa_series[iIm])/Sn;
			inv_ft[iIm] = (time_series[iRe]*spa_series[iIm] - time_series[iIm]*spa_series[iRe])/Sn;
		}
		gsl_fft_complex_radix2_backward(inv_ft, 1, data->NFFT);
		for (i=0; i<data->NFFT; i++)
		{
			if (fabs(inv_ft[2*i]) > max_corr)
			{
				max_corr = fabs(inv_ft[2*i]);
				tc_max = data->dt*(double)(i);
			} 
		}
	
		// recalculate SPA and its SNR
		eb->tc = -tc_max;
		snr_spa = 0.;
		for (i=1; i<=data->NFFT/2; i++)
		{
			if (i%data->under_samp == 0)
			{
				f = (double)(i)/(double)data->NFFT/data->dt;
				get_eccSPA(eb, f, &spaRe, &spaIm, j_max, j_min);
			
				iRe = 2*i;
				iIm = 2*i+1;

				spa_series[iRe] = spaRe;
				spa_series[iIm] = spaIm;
				
				Sn = get_Sn(f);
				
				snr_spa += (spa_series[iRe]*spa_series[iRe] + spa_series[iIm]*spa_series[iIm])/Sn;
			}
		}
		snr_spa *= 4./(double)data->NFFT/data->dt*(double)data->under_samp;
		snr_spa  = sqrt(snr_spa); 

		match = 0.;
		for (i=1; i<=data->NFFT/2; i++)
		{
			if (i%data->under_samp == 0)
			{
				iRe = 2*i;
				iIm = 2*i+1;
				f = (double)(i)/(double)data->NFFT/data->dt;
		
				Sn = get_Sn(f);

				match += (spa_series[iRe]*time_series[iRe] + spa_series[iIm]*time_series[iIm])/Sn;
				
			}
		}
		match *= 4./(double)data->NFFT*(double)data->under_samp/data->dt;
		match /= snr*snr_spa;
		
		eb->tc -= data->dt*(double)data->under_samp;
		
		for (m=0; m<m_max; m++)
		{
			for (i=1; i<=data->NFFT/2; i++)
			{
				if (i%data->under_samp == 0)
				{
					f = (double)(i)/(double)data->NFFT/data->dt;
					get_eccSPA(eb, f, &spaRe, &spaIm, j_max, j_min);
			
					iRe = 2*i;
					iIm = 2*i+1;

					spa_series[iRe] = spaRe;
					spa_series[iIm] = spaIm;
				}
			}
			snr_spa = get_SNR(spa_series, data);


			match = get_overlap(spa_series, time_series, data)/snr/snr_spa;
			eb->tc += 2.*data->dt/(double)m_max*(double)data->under_samp;
			l++;
			
			
			percent = (double)l/(double)((k_max)*(m_max+1));
			printProgress(percent);
			if (fabs(match)>max_match)
			{
				max_match = fabs(match);
				tc_mm = eb->tc - 2.*data->dt/(double)m_max;
				lc_mm = eb->lc;
			} 
		}

		eb->lc += PI2/(double)k_max;	
	}
	
	/////////////////////////////
	
	eb->lc = lc_mm;
	eb->tc = tc_mm;
	for (i=1; i<=data->NFFT/2; i++)
	{
		if (i%data->under_samp == 0)
		{
			f = (double)(i)/(double)data->NFFT/data->dt;
			get_eccSPA(eb, f, &spaRe, &spaIm, j_max, j_min);
			
			iRe = 2*i;
			iIm = 2*i+1;

			spa_series[iRe] = spaRe;
			spa_series[iIm] = spaIm;
		}
	}
	print_spa(spa_series, fopen("spa.dat", "w"), data);
	
	/////////////////////////////
	
	snr_spa = get_SNR(spa_series, data);
	fprintf(stdout, "\nFF: %f\n",   max_match);
	fprintf(stdout, "max tc: %e\n", tc_mm);
	fprintf(stdout, "max lc: %e\n", lc_mm);
	fprintf(stdout, "spa SNR: %f\n", snr_spa);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nFF duration: %f sec\n", time_spent);
	
	fprintf(stdout, "\n==============================================================\n");
	
	free(time_series);
	free(spa_series);
	free(inv_ft);
	
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
	
	eb->eLSO = 0.01;   // LSO eccentricity
	
	eb->F0 = 3.;

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
	
	data->under_samp = 100;
	data->NFFT     = (long)(pow(2, floor( log((double)(data->N))/log(2.)  ) + 1.));
	data->left_pad = (long)((data->NFFT - data->N)/2);
	
	fprintf(stdout, "Evolution and print runtime: %f sec\n", time_spent);
	fprintf(stdout, "\nN: %ld samples\n", data->N);
	fprintf(stdout, "FFT df: %lf\n", 1./data->NFFT/data->dt);
	fprintf(stdout, "TRU df: %lf\n", 1./data->N/data->dt);
	fprintf(stdout, "NFFT: %ld\n", data->NFFT);

	return;
}







