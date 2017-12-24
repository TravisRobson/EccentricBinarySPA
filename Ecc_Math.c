#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Ecc_SPA.h"
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