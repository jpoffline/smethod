
// fourier.cpp

#include "fourier.h"


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to compute Fourier transform

void ComputeFT(int n, double *datainput, fftw_complex *datainput_FT){

	// This routine computes the Fourier transform of REAL data
	// Input: datainput[] -- dim n, data to be Fourier transformed
	// Output datainput_FT[] -- dim n, contains UNnormlised Fourier mode coefficients
	
	
	// Declare plan to compute FT(datainput)
	fftw_plan p; 

	// Construct plan
	p = fftw_plan_dft_r2c_1d(n, datainput, datainput_FT, FFTW_ESTIMATE);	
	
	// Execute plan
	fftw_execute(p); 
	
	// Delete plan
	fftw_destroy_plan(p);
	
	
} // END ComputeFT()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to compute inverse Fourier transform

void ComputeiFT(int n, fftw_complex *datainput, double *datainput_iFT){

	// Declare plan to compute iFT(datainput)
	fftw_plan q;

	// Construct plan
	q = fftw_plan_dft_c2r_1d(n, datainput, datainput_iFT, FFTW_ESTIMATE);
	
	// Execute plan
	fftw_execute(q);
	
	// Delete plan
	fftw_destroy_plan(q);
	
} // END ComputeiFT()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF