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
void fill_spa_series(double *spa_series, struct EccBinary *eb, struct Data *data);
void fill_num_series(double *num_series, struct Data *data);
void fill_spa_series_new(double *spa_series, struct Data *data, double *spa_0, struct EccBinary *eb);
void fill_SPA_matrix(double **spa_matrix, struct EccBinary *eb, struct Data *data);
void spa_matrix_to_array(double **spa_matrix, double *spa_series, struct EccBinary *eb, struct Data *data);
void fill_Fisher(struct EccBinary *eb, struct Data *data);
void map_array_to_params(struct EccBinary *eb);
double invert_matrix(double **matrix, int N);

double find_max_tc(double *a, double *b, double *inv_ft, struct Data *data);


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

void map_array_to_params(struct EccBinary *eb)
{
	eb->Mc    = exp(eb->params[0])*TSUN;
	eb->F0    = exp(eb->params[1])*10.;
	eb->c0    = exp(eb->params[2])*1.;
	eb->lc    = eb->params[3];
	eb->tc    = -exp(eb->params[4])*10.;
	eb->R     = exp(eb->params[5])*(1.0e6*PC/C);
	eb->beta  = eb->params[6];
	eb->iota  = acos(eb->params[7]);
	eb->phi   = eb->params[8];
	eb->theta = acos(eb->params[9]);
	eb->psi   = eb->params[10];
	
	eb->c2beta = cos(2.*eb->beta);
	eb->s2beta = sin(2.*eb->beta);
	eb->ciota  = cos(eb->iota);
	eb->siota  = sin(eb->iota);
		
	return;
}

void fill_Fisher(struct EccBinary *eb, struct Data *data)
{

	return;
}

void spa_matrix_to_array(double **spa_matrix, double *spa_series, struct EccBinary *eb, struct Data *data)
{
	long i, j, iRe, iIm;
		
	double arg;
	
	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*(i*data->under_samp);
		iIm = 2*(i*data->under_samp)+1;

		spa_series[iRe] = 0.;
		spa_series[iIm] = 0.;
	}
		
	if (eb->lc == 0.)
	{
		for (j=0; j<eb->j_max; j++)
		{	
			for (i=1; i<=data->NFFT/2/data->under_samp; i++)
			{
				iRe = 2*(i*data->under_samp);
				iIm = 2*(i*data->under_samp)+1;
		
				spa_series[iRe] += spa_matrix[j][iRe];
				spa_series[iIm] += spa_matrix[j][iIm];
			}
		}
	}else 
	{
		for (j=0; j<eb->j_max; j++)
		{	
			arg = (double)(j+1)*eb->lc;
			for (i=1; i<=data->NFFT/2/data->under_samp; i++)
			{
				iRe = 2*(i*data->under_samp);
				iIm = 2*(i*data->under_samp)+1;
		
				spa_series[iRe] += spa_matrix[j][iRe]*cos(arg) - spa_matrix[j][iIm]*sin(arg);
				spa_series[iIm] += spa_matrix[j][iRe]*sin(arg) + spa_matrix[j][iIm]*cos(arg);
			}
		}
	}
	
	return;
}

void fill_SPA_matrix(double **spa_matrix, struct EccBinary *eb, struct Data *data)
{
	int temp_j_min, temp_j_max;
	
	long i, j, iRe, iIm;
	
	double f;
	double df = 1./(double)data->NFFT/data->dt;
	double spaRe, spaIm;
	
	temp_j_min = eb->j_min;
	temp_j_max = eb->j_max;
	
	for (j=0; j<temp_j_max; j++)
	{ 
		eb->j_min = j+1;
		eb->j_max = j+1;
		for (i=1; i<=data->NFFT/2/data->under_samp; i++)
		{
			iRe = 2*(i*data->under_samp);
			iIm = 2*(i*data->under_samp)+1;
		
			f = (double)(i*data->under_samp)*df;
		
			get_eccSPA(eb, f, &spaRe, &spaIm);

			spa_matrix[j][iRe] += spaRe;
			spa_matrix[j][iIm] += spaIm;
		}
	}
	eb->j_min = temp_j_min;
	eb->j_max = temp_j_max;

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
	eb->R  = 4.*1.0e6*PC/C; // 410 Mpc
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
	
	eb->j_max = 20;
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

void fill_spa_series(double *spa_series, struct EccBinary *eb, struct Data *data)
{
 	long i, iRe, iIm;
	
	double f;
	double df = 1./(double)data->NFFT/data->dt;
	double spaRe, spaIm;
	
	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*(i*data->under_samp);
		iIm = 2*(i*data->under_samp)+1;
		
		f = (double)(i*data->under_samp)*df;
		
		get_eccSPA(eb, f, &spaRe, &spaIm);
//fprintf(stdout, "test: %e\n", spaRe);
		spa_series[iRe] = spaRe;
		spa_series[iIm] = spaIm;
	}
	
	return;
}

double find_max_tc(double *a, double *b, double *inv_ft, struct Data *data)
{
	long i, iRe, iIm;
	
	double max_corr, tc_max, Sn, f;
	double df = 1./(double)data->NFFT/data->dt;
	
	max_corr = 0.;
	tc_max   = 0.;
	
	for (i=0; i<2*data->NFFT; i++)
	{
		inv_ft[i] = 0.;
	}
	for (i=1; i<=data->NFFT/2; i++)
	{
		f = (double)(i)*df;
		Sn = get_Sn(f);

		iRe = 2*i;
		iIm = 2*i+1;

		inv_ft[iRe] = (a[iRe]*b[iRe] + a[iIm]*b[iIm])/Sn;
		inv_ft[iIm] = (a[iRe]*b[iIm] - a[iIm]*b[iRe])/Sn;
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
	
	return tc_max;
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

void fill_spa_series_new(double *spa_series, struct Data *data, double *spa_0, struct EccBinary *eb)
{
	long i, iRe, iIm;
	
	double f;
	double df = 1./(double)data->NFFT/data->dt;
	double arg = PI2*eb->tc;
	
	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*i*data->under_samp;
		iIm = 2*(i*data->under_samp)+1;
		
		f = (double)(i*data->under_samp)*df;

		spa_series[iRe] =  spa_0[iRe]*cos(arg*f) + spa_0[iIm]*sin(arg*f);
		spa_series[iIm] = -spa_0[iRe]*sin(arg*f) + spa_0[iIm]*cos(arg*f);
	}

	return;
}

double invert_matrix(double **matrix, int N)
{
	int i,j;
	double cond;

	// Don't let errors kill the program (yikes)
	gsl_set_error_handler_off ();
	int err=0;

	// Find eigenvectors and eigenvalues
	gsl_matrix *GSLmatrix = gsl_matrix_alloc(N,N);
	gsl_matrix *GSLinvrse = gsl_matrix_alloc(N,N);
	gsl_matrix *cpy		  = gsl_matrix_alloc(N,N);
	gsl_matrix *SVDinv	  = gsl_matrix_alloc(N,N);
	gsl_matrix *Dmat	  = gsl_matrix_alloc(N,N);
	gsl_matrix *temp      = gsl_matrix_alloc(N, N);

	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			if(matrix[i][j]!=matrix[i][j])
			{
				fprintf(stdout, "error for parameters (%d, %d)\n", i, j);
				fprintf(stderr,"WARNING: nan matrix element, now what?\n");
			}
			gsl_matrix_set(GSLmatrix,i,j,matrix[i][j]);
			gsl_matrix_set(cpy,i,j,matrix[i][j]);
		}
	}

	//////
	//
	//	Calculate the SVD and condition number
	//
	///////

	gsl_matrix *V = gsl_matrix_alloc (N,N);
	gsl_vector *D = gsl_vector_alloc (N);
	gsl_vector *work = gsl_vector_alloc (N);

	gsl_linalg_SV_decomp(cpy, V, D, work);


	double max, min;
	max = -0.1;
	min = INFINITY;

	for (i=0; i<N; i++)
	{
		if (gsl_vector_get(D,i) > max) max = gsl_vector_get(D,i);

		if (gsl_vector_get(D,i) < min) min = gsl_vector_get(D,i);
	}

	cond = log10(max/min);
	
	
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++) 
		{
			if (i == j)
			{
				if (gsl_vector_get(D,i) < 1.0e-6) 
				{
					fprintf(stdout, "Near Singular value[%d]!!! ---> %e\n", i, gsl_vector_get(D,i));
					gsl_matrix_set(Dmat, i, j, 0.);
				} else
				{
					gsl_matrix_set(Dmat, i, j, 1./gsl_vector_get(D,i));
				}
				
			} else 
			{
				gsl_matrix_set(Dmat, i, j, 0.);
			}
		}
	
	}

	gsl_matrix_transpose(cpy);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Dmat, cpy,   0.0, temp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, temp, 0.0, SVDinv);
	
	
	////////

	gsl_permutation * permutation = gsl_permutation_alloc(N);

	err += gsl_linalg_LU_decomp(GSLmatrix, permutation, &i);
	err += gsl_linalg_LU_invert(GSLmatrix, permutation, GSLinvrse);

	if(err>0)
	{
		fprintf(stderr,"WARNING: singluar matrix\n");
		fflush(stderr);
	}else
	{
		//copy covariance matrix back into Fisher
		for(i=0; i<N; i++)
		{
			for(j=0; j<N; j++) 
			{
				matrix[i][j] = gsl_matrix_get(SVDinv, i, j);
			}
		}
	}

	gsl_vector_free(D);
	gsl_vector_free(work);
	gsl_matrix_free(V);
	gsl_matrix_free(Dmat);
	gsl_matrix_free(SVDinv);
	gsl_matrix_free(temp);

	gsl_matrix_free (GSLmatrix);
	gsl_matrix_free (GSLinvrse);
	gsl_permutation_free (permutation);

	return cond;
}



