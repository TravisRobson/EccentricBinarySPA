#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "Ecc_Math.h"
//#include "Ecc_SPA.h"
#include "Ecc_Binary.h"
#include "Constants.h"
#include "Detector.h"
#include "Ecc_IO.h"

double get_SNR(double *fft, struct Data *data)
{
	long i, iRe, iIm;
	
	double Sn, f;
	double snr = 0.;
	double df = 1./(double)data->NFFT/data->dt;

	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*i*data->under_samp;
		iIm = 2*(i*data->under_samp)+1;
		
		f = (double)(i*data->under_samp)*df;
		
		Sn = get_Sn(f);
		
		snr += (fft[iRe]*fft[iRe] + fft[iIm]*fft[iIm])/Sn;
	}

	snr *= 4.*df*(double)data->under_samp; 
	snr  = sqrt(snr); 	
	
	return snr;
}

double get_overlap(double *a, double *b, struct Data *data)
{
	long i, iRe, iIm;
	
	double Sn, f;
	double overlap = 0.;
	double df = 1./(double)data->NFFT/data->dt;
	
	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*i*data->under_samp;
		iIm = 2*(i*data->under_samp)+1;
		
		f = (double)(i*data->under_samp)*df;

		Sn = get_Sn(f);

		overlap += (a[iRe]*b[iRe] + a[iIm]*b[iIm])/Sn;
	}
	
	overlap *= 4.*df*(double)data->under_samp;
	
	return overlap;
}

double get_logL(double *a, double *b, struct Data *data)
{
	long i, iRe, iIm;
	
	double logL;
	
	double *dif;
	
	dif = malloc(2*data->NFFT*sizeof(double));
	
// 	for (i=0; i<2*data->NFFT; i++)
// 	{
// 		dif[i] = 0.;
// 	} 
	
	for (i=0; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*i*data->under_samp;
		iIm = 2*(i*data->under_samp)+1;
		
		dif[iRe] = a[iRe] - b[iRe];
		dif[iIm] = a[iIm] - b[iIm];
	}
	
	logL = -0.5*get_overlap(dif, dif, data);
	
	free(dif);

	return logL;
//	return 1.;
}

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N)
{
	int i,j;

	// Don't let errors kill the program (yikes)
	gsl_set_error_handler_off ();
	int err=0;

	// Find eigenvectors and eigenvalues
	gsl_matrix *GSLfisher = gsl_matrix_alloc(N,N);
	gsl_matrix *GSLcovari = gsl_matrix_alloc(N,N);
	gsl_matrix *GSLevectr = gsl_matrix_alloc(N,N);
	gsl_vector *GSLevalue = gsl_vector_alloc(N);

	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			if(matrix[i][j]!= matrix[i][j])fprintf(stderr,"WARNING: nan matrix element, now what?\n");
			gsl_matrix_set(GSLfisher, i, j, matrix[i][j]);
		}
	}

	// sort and put them into evec
	gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
	gsl_permutation * permutation = gsl_permutation_alloc(N);
	err += gsl_eigen_symmv (GSLfisher, GSLevalue, GSLevectr, workspace);
	err += gsl_eigen_symmv_sort (GSLevalue, GSLevectr, GSL_EIGEN_SORT_ABS_ASC);

	// eigenvalues destroy matrix
	for(i=0; i<N; i++) for(j=0; j<N; j++) gsl_matrix_set(GSLfisher, i, j, matrix[i][j]);

	err += gsl_linalg_LU_decomp(GSLfisher, permutation, &i);
	err += gsl_linalg_LU_invert(GSLfisher, permutation, GSLcovari);

	if(err>0)
	{
		for(i=0; i<N; i++)for(j=0; j<N; j++)
		{
			evector[i][j] = 0.0;
			if(i==j)
			{
				evector[i][j]=1.0;
				evalue[i]=1./matrix[i][j];
			}
		}

	}
	else
	{
		//unpack arrays from gsl inversion
		for(i=0; i<N; i++)
		{
			evalue[i] = gsl_vector_get(GSLevalue, i);
			for(j=0; j<N; j++)
			{
				evector[i][j] = gsl_matrix_get(GSLevectr, i, j);
				if(evector[i][j] != evector[i][j]) evector[i][j] = 0.;
			}
		}
		
		//copy covariance matrix back into Fisher
		for(i=0; i<N; i++)
		{
			for(j=0; j<N; j++)
			{
				matrix[i][j] = gsl_matrix_get(GSLcovari, i, j);
			}
		}

		//cap minimum size eigenvalues
		for(i=0; i<N; i++)
		{
			if(evalue[i] != evalue[i] || evalue[i] <= 10.0) evalue[i] = 10.;
			//fprintf(stdout, "here\n");
		}
	}

	gsl_vector_free (GSLevalue);
	gsl_matrix_free (GSLfisher);
	gsl_matrix_free (GSLcovari);
	gsl_matrix_free (GSLevectr);
	gsl_eigen_symmv_free (workspace);
	gsl_permutation_free (permutation);
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

double max_spa_tc_lc(struct EccBinary *eb, struct Data *data, double *spa_series, double *num_series, double snr)
{
	int k_max = 150; // for lc search (TODO: SEEM TO NEED LESS FOR LESS ECCENTRIC SYSTEMS)
	int l     = 0;
	
	long i, j, k;
	
	double match, snr_spa, percent;
	double max_match = 0.;
	double tc_mm     = 0.;
	double lc_mm     = 0.;
	
	double *spa_0, *inv_ft;
	double **spa_mat;	
		
	spa_0   = malloc(2*data->NFFT*sizeof(double));
	inv_ft  = malloc(2*data->NFFT*sizeof(double));
	spa_mat = malloc(eb->j_max*sizeof(double *));
	for (j=0; j<eb->j_max; j++) spa_mat[j] = malloc(2*data->NFFT*sizeof(double));
	
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_0[i]  = 0.;
		for (j=0; j<eb->j_max; j++) spa_mat[j][i] = 0.;
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

	// free memory
	
	free(spa_0);
	free(inv_ft);
	for (j=0; j<eb->j_max; j++) free(spa_mat[j]);

	return max_match;
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

void fill_spa_series(double *spa_series, struct EccBinary *eb, struct Data *data)
{
 	long i, iRe, iIm;
	
	double f;
	double df = 1./(double)data->NFFT/data->dt;
	double spaRe, spaIm;
	
//	struct EccBinary eb_local = *eb;
	
	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*(i*data->under_samp);
		iIm = 2*(i*data->under_samp)+1;
		
		f = (double)(i*data->under_samp)*df;
		
//		get_eccSPA(&eb_local, f, &spaRe, &spaIm);
		get_eccSPA(eb, f, &spaRe, &spaIm);

		spa_series[iRe] = spaRe;
		spa_series[iIm] = spaIm;
	}
	
	return;
}












