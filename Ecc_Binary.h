#ifndef Ecc_Binary_h
#define Ecc_Binary_h

#include <gsl/gsl_spline.h>

struct EccBinary
{
	gsl_interp_accel *acc;
	gsl_spline *spline;
	
	gsl_interp_accel *e_acc;
	gsl_spline *e_spline;
	
	double *params;
	int NP;
	double **Fisher;
	
	double Jj, Jjm1, Jjp1;
	double inve;
	
	int j_min, j_max;
	
	double m1, m2;		    // component masses of binary
	double m, eta, Mc, mu;  // total mass, symmetric mass ratio, chirp mass, reduced mass
	
	double iota, beta;	    // angles associated with orbit orientation
	double theta, phi;		// angles associated with sky location
	double Fp, Fc;		    // antenna patterns for plus (p) and cross (c) polarizations
	double psi;				// polarization angle
	
	// convenient quantities for SPA
	double c2beta, s2beta;  
	double ciota, siota;
	
	double c0;				// constant of integration for Forb, e equations
	double eLSO, FLSO;      // last stable orbit (LSO) eccentricity and mean orbital frequency
	double rhoLSO;			// dimensional semi-major axis as LSO
	
	double R;				// distance to source
	
	double tc, lc;			// constants of integration in phase of the system
	
	// for the numerical evolution part
	double F0, e0, p0;
};

void set_m(struct EccBinary *eb);
void set_eta(struct EccBinary *eb);
void set_mu(struct EccBinary *eb);
void set_Mc(struct EccBinary *eb);

void set_Antenna(struct EccBinary *eb);

double calc_sigma(double e);

void set_FLSO(struct EccBinary *eb);
void set_c0(struct EccBinary *eb);
double solve_sigma_e(struct EccBinary *eb);

double get_e(double x);

double get_hp(struct EccBinary *eb, double phi, double p, double e);
double get_hc(struct EccBinary *eb, double phi, double p, double e);


#endif /* Ecc_Binary_h */