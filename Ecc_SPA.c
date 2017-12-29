#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_result.h>

#include "Ecc_SPA.h"
#include "Ecc_Binary.h"
#include "Constants.h"
#include "Ecc_IO.h"
#include "Ecc_Math.h"

double get_Cpj(struct EccBinary *eb, double j, double e)
{
	double Jj, Jjm1, Jjp1;
	double a, b, c;
	double c2b, ci, si; 
	
	c2b = eb->c2beta;
	ci  = eb->ciota;
	si  = eb->siota;

	Jj   = eb->Jj;
	Jjm1 = eb->Jjm1;
	Jjp1 = eb->Jjp1;
	
	a = c2b*(1. + ci*ci)*e*(e*e - 1.)*j;
	b = c2b*(1. + ci*ci)*(e*e - 2.) + e*e*si*si;
	c = c2b*(1. + ci*ci)*e*(1. - e*e)*j;
		
	return 2.*(a*Jjm1 - b*Jj + c*Jjp1)*eb->inve*eb->inve;
}

double get_Spj(struct EccBinary *eb, double j, double e)
{
	double Jj, Jjm1;
	double a, b;
	double ci, s2b;
	
	ci  = eb->ciota;
	s2b = eb->s2beta;

	Jj   = eb->Jj;
	Jjm1 = eb->Jjm1;
	
	a = e;
	b = (1. + (1. - e*e)*j);
	
	return 4.*(1. + ci*ci)*s2b*sqrt(1. - e*e)*(a*Jjm1 - b*Jj)*eb->inve*eb->inve;
}

double get_Ccj(struct EccBinary *eb, double j, double e)
{
	double Jj, Jjm1;
	double a, b;
	double s2b, ci;
	
	ci  = eb->ciota;
	s2b = eb->s2beta;

	Jj   = eb->Jj;
	Jjm1 = eb->Jjm1;

	a = 2.*e*(1. - e*e)*j;
	b = e*e - 2. + 2.*(e*e - 1.)*j;
	
	return 4*s2b*ci*(a*Jjm1 + b*Jj)*eb->inve*eb->inve;
}

double get_Scj(struct EccBinary *eb, double j, double e)
{
	double Jj, Jjm1;
	double a, b;
	double c2b, ci;
	
	c2b = eb->c2beta;
	ci  = eb->ciota;

	Jj   = eb->Jj;
	Jjm1 = eb->Jjm1;
	
	a = e;
	b = -1. + (e*e - 1.)*j;
	
	return 8.*c2b*ci*sqrt(1. - e*e)*(a*Jjm1 + b*Jj)*eb->inve*eb->inve;
}

double get_es(struct EccBinary *eb, double alpha)
{
	// bisection method implementation
	double RHS;			// RHS of equation to be inverted
	double error, goal; // error, goal precision
		
	long itr_max = 100000;
	long itr = 0;
	
	goal = 1.0e-10; // error goal, absolute difference
	
	error = 1.;		// absolute error
	
	RHS    = pow(alpha, 2./3.);
	
	double a, b;
	double c1, fc1;
	
	a = 0.;
	b = 0.99999999;
	
	while (error > goal)
	{	
		c1  = 0.5*(a + b);
		
		fc1 = calc_sigma(c1) - RHS;
		error = fabs( calc_sigma(c1) - RHS );
		if (fc1 < 0.) a = c1;
		else b = c1;
		
		itr += 1;
		if (itr >= itr_max) 
		{	// ensure that we don't get stuck in a infinite loop
			fprintf(stdout, "Max iterations exceeded in 'get_es()' \n");
			break;
		}
	}

	return c1;
}

double interp_es(struct EccBinary *eb, double alpha)
{
	//return gsl_spline_eval(eb->e_spline, alpha, eb->e_acc);
	if (alpha < 0.000068257 || alpha > 13152.153520857013)
	{
		return get_es(eb, alpha);
	} else 
	{
		return gsl_spline_eval(eb->e_spline, alpha, eb->e_acc);
	}
	
}

double term(int r, double x)
{
	double result;
	
	result  = gsl_sf_gamma(1. + r)/gsl_sf_gamma(1.);
	result *= gsl_sf_gamma(-1181./2299 + r)/gsl_sf_gamma(-1181./2299.);
	result *= gsl_sf_gamma(3./2. + r)/gsl_sf_gamma(3./2.);
	result *= gsl_sf_gamma(43./19. - 1. + r)/gsl_sf_gamma(43./19. - 1.);
	result /= gsl_sf_gamma(43./19. + r - 1. + r)/gsl_sf_gamma(43./19. + r - 1.);
	result /= gsl_sf_gamma(43./19. + 2.*r)/gsl_sf_gamma(43./19.);
	result /= gsl_sf_fact(r);
	
	result *= pow((121.*x)/(304. + 121.*x), r)*pow(x/(x - 1.), r);


	result *= gsl_sf_hyperg_2F1(1. + (double)r, -1181./2299. + (double)r, 43./19. + 2.*(double)r, 121.*x/(304. + 121.*x)); 
	//result *= gsl_sf_hyperg_2F1(1. + (double)r, 3./2. + (double)r, 43./19. + 2.*(double)r, x/(x - 1.));
	result *= pow(1. - x, 1. + (double)r);
	result *= gsl_sf_hyperg_2F1(1. + (double)r, 43./19. - 3./2. + (double)r, 43./19. + 2.*(double)r, x);
	
	return result;
}

double get_It(struct EccBinary *eb, double e)
{
	int i, n;
	double x, result, temp;
	
	n = 10;
	
// 	     if (e < 0.1) n = 2;
// 	else if (e < 0.3) n = 3;
// 	else if (e < 0.4) n = 4;
// 	else if (e < 0.6) n = 5;
// 	else if (e < 0.7) n = 6;
// 	else if (e < 0.8) n = 7;
// 	else if (e < 0.9) n = 9;
	
	x = e*e;
	
	result  = 19./48.*pow(e, 48./19.)*pow(1. + 121./304.*x, 1181./2299.)/pow(1. - x, 3./2.);

	temp = 0.;
	
	for (i=0; i<n; i++)
	{
		temp += term(i, x);
	}	
		
	return result*temp;
}

double get_Il(struct EccBinary *eb, double e)
{
	//fprintf(stdout, "fuck it all: %e\n", e);
	return 19./30.*pow(e, 30./19.)*gsl_sf_hyperg_2F1(124./2299., 15./19., 34./19., -121./304.*e*e) ;
}

double get_phase(struct EccBinary *eb, int j, double e, double f)
{
	double result;
	
	result  = 15./304./pow(eb->c0, 8./3.)/pow(PI2*eb->Mc, 5./3.); // fprintf(stdout, "result: %e\n", eb->Mc);
	result *= f*get_It(eb, e) - (double)(j)*eb->c0*get_Il(eb, e);
	result += -PI2*f*eb->tc + (double)(j)*eb->lc;
	
	return result;
}	

double interp_phase(struct EccBinary *eb, double j, double alpha, double f)
{
	double result;
	
	result  = 15./304./pow(eb->c0, 8./3.)/pow(PI2*eb->Mc, 5./3.)*f;
	result *= gsl_spline_eval(eb->spline, alpha, eb->acc);
	result += -PI2*f*eb->tc + j*eb->lc;

	return result;
}

double get_Aj(struct EccBinary *eb, double j, double e)
{
	double result;
	
	result = pow(j, 2./3.)*pow(1. - e*e, 7./4.)/sqrt(1. + 73./24.*e*e + 37./96.*e*e*e*e);
	
	return result;
}

void get_harmonic(struct EccBinary *eb, int j, double f, double *hjRe, double *hjIm)
{
	double es;	   // stationary eccentricity
	double amp;    // amplitude
	double phase;  // phase
	double Cpj, Spj, Ccj, Scj;
	
	double jd = (double)j;
	double alpha;
	double e;
	double K1, K2, dJj;
	
	if (f > jd*eb->FLSO || f < jd*eb->F0)
	{
		*hjRe = 0.;
		*hjIm = 0.;
		return;
	} 

	alpha = eb->c0*jd/f;
	es = interp_es(eb, alpha);	
	eb->inve = 1./es;
	e = es;


	if (e < 0.8){
		eb->Jj   = gsl_sf_bessel_Jn(j,   jd*e);
		eb->Jjm1 = gsl_sf_bessel_Jn(j-1, jd*e);
		eb->Jjp1 = gsl_sf_bessel_Jn(j+1, jd*e);
	}else {
	
		K1 = gsl_sf_bessel_Kn(1./3., jd*(-sqrt(1 - pow(e, 2.)) + log((1. + sqrt(1. - pow(e, 2.)))*eb->inve)));
		K2 = gsl_sf_bessel_Kn(2./3., jd*(-sqrt(1 - pow(e, 2.)) + log((1. + sqrt(1. - pow(e, 2.)))*eb->inve)));	
		
		eb->Jj = pow(2., 1./3.)*pow(3.,1./6.)*pow(pow(-sqrt(1 - pow(e, 2.)) 
				 + log((1. + sqrt(1. - pow(e, 2.)))/e), 2./3.)
	         	 /(1 - pow(e,2)),0.25)*((K1*pow(-sqrt(1 - pow(e,2)) 
	         	 + log((1. + sqrt(1 - pow(e,2)))/e),1./3.))/(pow(2,1./3.)*pow(3,1./6.)*M_PI) 
	         	 + (K2*pow(-sqrt(1 - pow(e,2)) + log((1 + sqrt(1 - pow(e,2)))/e),1./3.)*((-5/(3.*pow(1 - pow(e,2),1.5)) 
	         	 + 1./sqrt(1 - pow(e,2)))/8. + 5/(72.*(-sqrt(1 - pow(e,2)) 
	         	 + log((1 + sqrt(1 - pow(e,2)))/e)))))/(pow(2,1./3.)*pow(3,1./6.)*jd*M_PI));
		dJj    = -(((1 - pow(e,2))*(-gsl_sf_bessel_Kn(2./3.,j*(-sqrt(1 - pow(e,2)) 
	         	 + log((1 + sqrt(1 - pow(e,2)))/e))) 
	         	 - (gsl_sf_bessel_Kn(1./3.,j*(-sqrt(1 - pow(e,2)) + log((1 + sqrt(1 - pow(e,2)))/e)))
	         	 *((-6 + 27*pow(e,2))/pow(1 - pow(e,2),1.5) + 7/(sqrt(1 - pow(e,2)) 
	         	 - log((1 + sqrt(1 - pow(e,2)))/e))))/(72.*j)))/(e*M_PI*pow((1 - pow(e,2))/pow(-sqrt(1 - pow(e,2)) 
	         	 + log((1 + sqrt(1 - pow(e,2)))/e),2./3.),0.75)));     
	           
	eb->Jjm1 = eb->Jj/e + dJj;
	           
	eb->Jjp1 = eb->Jj/e - dJj;
	} 
	
	amp  = get_Aj(eb, jd, es);
	
	//phase =  get_phase(eb, j, es, f); 	
	phase = interp_phase(eb, j, alpha, f);
	
	Cpj = get_Cpj(eb, jd, es);
	Spj = get_Spj(eb, jd, es);
	Ccj = get_Ccj(eb, jd, es);
	Scj = get_Scj(eb, jd, es);
	
	*hjRe = amp*( (eb->Fp*Cpj + eb->Fc*Ccj)*cos(phase + PIon4) 
	             +(eb->Fp*Spj + eb->Fc*Scj)*sin(phase + PIon4) );

	*hjIm = amp*(-(eb->Fp*Spj + eb->Fc*Scj)*cos(phase + PIon4) 
	             +(eb->Fp*Cpj + eb->Fc*Ccj)*sin(phase + PIon4) );	

	return;
}

// calculate the real an imaginary parts of the SPA at a given f
void get_eccSPA(struct EccBinary *eb, double f, double *spaRe, double *spaIm)
{
	int j;
	int j_min = eb->j_min;
	int j_max = eb->j_max;
	
	double re, im;
	double a;
	
	*spaRe = 0.;
	*spaIm = 0.;
	
	for (j=j_min; j<=j_max; j++)
	{
		get_harmonic(eb, j, f, &re, &im);
		
		*spaRe += re;
		*spaIm += im;
	}
	
	a = 0.5*sqrt(5.*M_PI/48.)*pow(eb->Mc, 5./6.)/eb->R*pow(PI2*f, -7./6.);
	*spaRe *= a;
	*spaIm *= a;

	return;
}

void setup_interp(struct EccBinary *eb)
{
	int i;
	int n = 10000;
	
	double *x, *psi, *e_vec;
	
	FILE *in_file;
	
	gsl_interp_accel *acc      = gsl_interp_accel_alloc();
	gsl_spline       *spline   = gsl_spline_alloc(gsl_interp_cspline, n);
	gsl_interp_accel *e_acc    = gsl_interp_accel_alloc();
	gsl_spline       *e_spline = gsl_spline_alloc(gsl_interp_cspline, n);
	
	x     = malloc(n*sizeof(double));
	psi   = malloc(n*sizeof(double));
	e_vec = malloc(n*sizeof(double));
	
	// Interpolate Psi(alpha)
	for(i=0; i<n; i++)
	{
		x[i]   = 0.;
		psi[i] = 0.;
	}	
	i = 0;
	in_file = fopen("psi_tbl.dat", "r");
	while ( !feof(in_file) )
	{	
		fscanf(in_file, "%lf%lf\n", &x[i], &psi[i]);
		i ++;
	}
	fclose(in_file);
	
	eb->spline = spline;
	eb->acc    = acc;
	gsl_spline_init(eb->spline, x, psi, n);

	
	// Interpolate e(alpha)
	for(i=0; i<n; i++)
	{
		x[i]   = 0.;
		e_vec[i] = 0.;
	}	
	i = 0;
	in_file = fopen("e_tbl.dat", "r");
	while ( !feof(in_file) )
	{	
	
		fscanf(in_file, "%lf%lf\n", &x[i], &e_vec[i]);
		i += 1;
	}
	fclose(in_file);

	eb->e_spline = e_spline;
	eb->e_acc = e_acc;
	gsl_spline_init(eb->e_spline, x, e_vec, n);
		
	free(x);
	free(psi);
	free(e_vec);	
	
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
	
	for (i=1; i<=data->NFFT/2/data->under_samp; i++)
	{
		iRe = 2*(i*data->under_samp);
		iIm = 2*(i*data->under_samp)+1;
		
		f = (double)(i*data->under_samp)*df;
		
		get_eccSPA(eb, f, &spaRe, &spaIm);

		spa_series[iRe] = spaRe;
		spa_series[iIm] = spaIm;
	}
	
	return;
}

void calc_Fisher(struct EccBinary *eb, struct Data *data)
{	
	long iRe, iIm; 
	long i, k;
	
	double ep = 1.0e-6;
	
	double *spa_p, *spa_m, **spa_dif;
	
	struct EccBinary *eb_p = malloc(sizeof(struct EccBinary));
	struct EccBinary *eb_m = malloc(sizeof(struct EccBinary));
	
	eb_p->NP = eb->NP;
	eb_m->NP = eb->NP;
	
	eb->Fisher = malloc(eb->NP*sizeof(double *));
	for (i=0; i<eb->NP; i++) eb->Fisher[i] = malloc(eb->NP*sizeof(double));
	
	eb_p->params = malloc(eb_p->NP*sizeof(double));
	eb_m->params = malloc(eb_m->NP*sizeof(double));
	
	spa_p   = malloc(2*data->NFFT*sizeof(double));
	spa_m   = malloc(2*data->NFFT*sizeof(double));
	spa_dif = malloc(eb->NP*sizeof(double *));
	for (i=0; i<eb->NP; i++) spa_dif[i] = malloc(2*data->NFFT*sizeof(double));

	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<eb->NP; k++) eb->Fisher[i][k] = 0.;
	}
	for (i=0; i<eb->NP; i++)
	{
		eb_p->params[i] = 0.;
		eb_m->params[i] = 0.;
	}
	for (i=0; i<2*data->NFFT; i++)
	{
		spa_p[i]   = 0.;
		spa_m[i]   = 0.;
	}
	for (i=0; i<eb->NP; i++)
	{
		for (k=0; k<2*data->NFFT; k++) spa_dif[i][k] = 0.;
	}

	setup_interp(eb_p);
	setup_interp(eb_m);
	
	eb_p->FLSO  = eb->FLSO;
	eb_m->FLSO  = eb->FLSO;
	eb_p->j_min = eb->j_min;
	eb_m->j_min = eb->j_min;
	eb_p->j_max = eb->j_max;
	eb_m->j_max = eb->j_max;

	for (i=0; i<eb->NP; i++)
	{
		eb_p->params[i] = eb->params[i];
		eb_m->params[i] = eb->params[i];
	}
	map_array_to_params(eb_p);
	map_array_to_params(eb_m);
	set_Antenna(eb_p);
	set_Antenna(eb_m);
	
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
	for (i=0; i<eb->NP; i++) free(spa_dif[i]);
	
	// use symmetry to populate lower triangle
	for (i=0; i<eb->NP; i++)
	{
		for (k=i; k<eb->NP; k++) eb->Fisher[k][i] = eb->Fisher[i][k];

	}		

	return;
}











