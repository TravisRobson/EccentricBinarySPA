#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Ecc_Binary.h"
#include "Constants.h"


void set_m(struct EccBinary *eb)
{
	eb->m = eb->m1 + eb->m2;
	return;
}

void set_eta(struct EccBinary *eb)
{
	eb->eta = eb->m1*eb->m2/(eb->m*eb->m);
	return;
}

void set_mu(struct EccBinary *eb)
{
	eb->mu = eb->m1*eb->m2/(eb->m);
	return;
}

void set_Mc(struct EccBinary *eb)
{
	eb->Mc = pow(eb->m1*eb->m2, 3./5.)/pow(eb->m, 1./5.);
	return;
}

void set_Antenna(struct EccBinary *eb)
{
	double ctheta, c2phi, c2psi, s2phi, s2psi;
	
	ctheta = cos(eb->theta);
	c2phi  = cos(2.*eb->phi);
	s2phi  = sin(2.*eb->phi);
	c2psi  = cos(2.*eb->psi);
	s2psi  = sin(2.*eb->psi);
	
	eb->Fp = 0.5*(1. + ctheta*ctheta)*c2phi*c2psi - ctheta*s2phi*s2psi;
	eb->Fc = 0.5*(1. + ctheta*ctheta)*c2phi*s2psi + ctheta*s2phi*c2psi;
	
	return;
}

double calc_sigma(double e)
{
	double sigma;
	
	sigma  = pow(e, 12./19.)/(1. - e*e);
	sigma *= pow(1. + 121./304.*e*e, 870./2299.);
	
	return sigma;
}	

void set_FLSO(struct EccBinary *eb)
{	
	eb->rhoLSO = (6. + 2.*eb->eLSO)/(1. + eb->eLSO);
	
	//eb->FLSO = pow(1. - eb->eLSO, 3./2.)/PI2/eb->m/pow(eb->rhoLSO, 3./2.);
	
	// From PC
	eb->FLSO = pow((1.+eb->eLSO)/(6. + 2.*eb->eLSO), 3./2.)/eb->m/PI2;

	return;
}

void set_c0(struct EccBinary *eb)
{
	eb->c0  = pow(calc_sigma(eb->eLSO), 3./2.)*pow(1. - eb->eLSO, 3./2.);
	eb->c0 /= pow(eb->rhoLSO, 3./2.)*PI2*eb->m;

	return;
}

double solve_sigma_e(struct EccBinary *eb)
{
	// bisection method implementation
	double RHS;			// RHS of equation to be inverted
	double error, goal; // error, goal precision
		
	long itr_max = 100000;
	long itr = 0;
	
	goal = 1.0e-12; // error goal, absolute difference
	
	error = 1.;		// absolute error
	
	RHS    = pow(eb->c0/eb->F0, 0.66666666666667);
	
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

double get_e(double x)
{
	// bisection method implementation
	double RHS;			// RHS of equation to be inverted
	double error, goal; // error, goal precision
		
	long itr_max = 100000;
	long itr = 0;
	
	goal = 1.0e-10; // error goal, absolute difference
	
	error = 1.;		// absolute error
	
	RHS    = pow(x, 0.66666666666667);
	
	double a, b;
	double c1, fc1;
	
	a = 0.;
	b = 0.9999999999;
	
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

double get_hp(struct EccBinary *eb, double phi, double p, double e)
{
	double mu, R;
	double ciota, siota, beta;
	double hp;
	
	mu    = eb->mu;
	R     = eb->R;
	ciota = eb->ciota;
	siota = eb->siota;
	beta  = eb->beta;
	
	hp  = 2.*cos(2.*(phi - beta)) +5.*e/2.*cos(phi - 2.*beta);
	hp += 0.5*e*cos(3.*phi - 2.*beta) + e*e*cos(2.*beta);
	hp *= (1. + ciota*ciota);
	hp += (e*cos(phi) + e*e)*siota*siota;
	hp *= -mu/p/R;

	return hp;
}


double get_hc(struct EccBinary *eb, double phi, double p, double e)
{
	double mu, R;
	double ciota, beta;
	double hc;
	
	mu    = eb->mu;
	R     = eb->R;
	ciota = eb->ciota;
	beta  = eb->beta;
	
	hc  = 4.*sin(2.*(phi - beta)) + 5.*e*sin(phi - 2.*beta);
	hc += e*sin(3.*phi - 2.*beta) - 2.*e*e*sin(2.*beta);
	hc *= ciota;
	hc *= -mu/p/R;
	
	return hc;
}


























