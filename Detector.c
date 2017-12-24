#include <stdio.h>
#include <math.h>
#include <stdlib.h>




double get_Sn(double f)
{
	// Taken from equation (2.1) in Cutler Flanagan '94
	double f0;
	double S0, Sn;
	
	double fonf0;
	
	f0 = 70.; 	 // (Hz)
	S0 = 3.e-48; // Hz^{-1}
	fonf0 = f*f/f0/f0;
	
	if (f < 10.) return INFINITY;
	else
	{
		Sn = 1./fonf0/fonf0 + 2.*(1. + fonf0);
		Sn *= S0;
	}
	
	return Sn;
}

