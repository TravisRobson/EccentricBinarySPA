#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

#include "Ecc_Binary.h"
#include "Constants.h"
#include "Ecc_IO.h"

int func(double t, const double y[], double f[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	
	struct EccBinary *eb = (struct EccBinary *)params;
	double m =  eb->m;
	double mu = eb->mu;


	// orbital phase evolution
	double what;
	what = pow(y[1], 3./2.);
	
	f[0]  = (1. + y[2]*cos(y[0]))*(1. + y[2]*cos(y[0]))/what/m;
	
	// semi-latus rectum evolution
	f[1]  = -64./5.*mu/m/m*pow(1. - y[2]*y[2], 3./2.)/y[1]/y[1]/y[1];
	f[1] *= (1. + 7./8.*y[2]*y[2]);
	
	// eccentricity evolution
	f[2]  = -304./15.*mu/m/m*pow(1. - y[2]*y[2], 3./2.)/y[1]/y[1]/y[1]/y[1];
	f[2] *= y[2]*(1. + 121./304*y[2]*y[2]);
	
	return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	double m00, m01, m02;
	double m10, m11, m12;
	double m20, m21, m22;
	
	struct EccBinary *eb = (struct EccBinary *)params;
	double m =  eb->m;
	double mu = eb->mu;
  
  	double what;
	
	what = pow(y[1], 3./2.);
  
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
	
	gsl_matrix * mat = &dfdy_mat.matrix; 
	

	m00  = -2.*y[2]*sin(y[0])*(1. + y[2]*cos(y[0]))/what/m;
	m01  = -3./2./y[1]*(1. + y[2]*cos(y[0]))*(1. + y[2]*cos(y[0]))/what/m;
	m02  = 2.*(1. + y[2]*cos(y[0]))*cos(y[0])/what/m;
	
	m10  = 0.;
	m11  = -3*64./5.*mu/m/m*pow(1. - y[2]*y[2], 3./2.)/y[1]/y[1]/y[1]/y[1]*(1. + 7./8.*y[2]*y[2]);
	m12  = -64./5.*mu/m/m*pow(1. - y[2]*y[2], 3./2.)/y[1]/y[1]/y[1];
	m12 *= -3.*y[2]/(1. - y[2]*y[2])*(1. + 7./8.*y[2]*y[2]) + 7./4.*y[2];
	
	m20  = 0.;
	m21  = 4.*304./15.*mu/m/m*pow(1. - y[2]*y[2], 3./2.)/y[1]/y[1]/y[1]/y[1];
	m21 *= y[2]*(1. + 121./304*y[2]*y[2])/y[1];
	m22  = -304./15.*mu/m/m*pow(1. - y[2]*y[2], 3./2.)/y[1]/y[1]/y[1]/y[1];
	m22 *= (1. + 121./304.*y[2]*y[2] + 121./152.*y[2]*y[2]);
	
	// first row
	gsl_matrix_set(mat, 0, 0, m00);
	gsl_matrix_set(mat, 0, 1, m01);
	gsl_matrix_set(mat, 0, 1, m02);		
	
	// second row
	gsl_matrix_set(mat, 1, 0, m10);
	gsl_matrix_set(mat, 1, 1, m11);
 	gsl_matrix_set(mat, 1, 2, m12);

 	// third row
	gsl_matrix_set(mat, 2, 0, m20);
	gsl_matrix_set(mat, 2, 1, m21);
	gsl_matrix_set(mat, 2, 2, m22);


	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	
	return GSL_SUCCESS;
}


void evolve_binary(struct EccBinary *eb, struct Data *data, double y[], gsl_odeiv2_driver *d, FILE *fptr)
{
	int status;
	
	long i;
	
	double t, ti, Forb;
	double dt, h;
	
	dt = data->dt;
	i = 0; t = 0.;
	
	while (Forb < eb->FLSO)
	{
		ti = (double)i*dt;
		status = gsl_odeiv2_driver_apply(d, &t, ti, y);
		
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		
		Forb = pow((1. - y[2]*y[2])/y[1], 3./2.)/PI2/eb->m;
		h = get_hp(eb, y[0], y[1], y[2])*eb->Fp + get_hc(eb, y[0], y[1], y[2])*eb->Fc;

		fprintf(fptr, "%e %e %e %e %e %e\n", t, y[0], y[1], y[2], Forb, h);
		i += 1;
	}
	fclose(fptr);
	
	data->N = i;	// number of data samples

	return;
}


void evolve_no_RR(struct EccBinary *eb, struct Data *data, double y[], gsl_odeiv2_driver *d, FILE *fptr)
{
	int status;
	
	long i;
	
	double t, ti;
	double dt, h;
	
	dt = data->dt;
	i = 0; t = 0.;
	
	while (ti < 40.)
	{
		ti = (double)i*dt;
		status = gsl_odeiv2_driver_apply(d, &t, ti, y);
		
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		
		h = get_hp(eb, y[0], eb->p0, eb->e0)*eb->Fp + get_hc(eb, y[0], eb->p0, eb->e0)*eb->Fc;

		fprintf(fptr, "%e %e %e\n", t, y[0], h);
		i += 1;
	}
	fclose(fptr);
	
	data->N = i;	// number of data samples

	return;
}

int func_no_RR(double t, const double y[], double f[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	
	struct EccBinary *eb = (struct EccBinary *)params;
	double m  = eb->m;
	double p0 = eb->p0;
	double e0 = eb->e0;
	

	// orbital phase evolution
	double what;
	what = pow(p0, 3./2.);
	
	f[0]  = (1. + e0*cos(y[0]))*(1. + e0*cos(y[0]))/what/m;
	
	return GSL_SUCCESS;
}

int jac_no_RR(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	double m00;
	
	struct EccBinary *eb = (struct EccBinary *)params;
	double m =  eb->m;
	double p0 = eb->p0;
	double e0 = eb->e0;
  
  	double what;
	
	what = pow(p0, 3./2.);
  
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
	gsl_matrix * mat = &dfdy_mat.matrix; 

	m00  = -2.*e0*sin(y[0])*(1. + e0*cos(y[0]))/what/m;
	gsl_matrix_set(mat, 0, 0, m00);

	dfdt[0] = 0.0;
	
	return GSL_SUCCESS;
}




















