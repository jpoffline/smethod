
// derivs_FFT.cpp

// This contains routines to compute all spatial derivatives
// via Fast Fourier Transform (FFT).

#include "derivs_FFT.h"

void ComputeLaplacian_FFT(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
	// Routine to compute Laplacian of "field" via FFT
	
	// Temporary field array
	double *tempfield = new double[params->imax];
	// Temporary array to hold laplacian of a given component
	double *FFTlap_tempfield = new double[params->imax];
	
	for(int com = 0; com < params->cmax; com++){
	
		// Extract the component of the field.
		// Dump into temp-field array.
		for(int i = 0; i < params->imax; i++){
			tempfield[i] = field->vals[ field->ind(grid->now,com,i,grid,field) ];
		}
		
		// Compute laplacian(tempfield) via FFT
		FFT_DDeriv(params->imax, params->h, tempfield, FFTlap_tempfield);
		
		// Dump result into object in field struct
		for(int i = 0; i < params->imax; i++){
			field->FFTlap[ field->ind(grid->now,com,i,grid,field) ] = FFTlap_tempfield[i];
		}
		
	} // END com-loop
	
	delete tempfield;
	delete FFTlap_tempfield;

} // END ComputeLaplacian()


void FFT_Deriv(int n, double h, double *datainput, double *datainput_derivFT){

	// Function to return the first derivative
	// of real data, computed via FFT
	
	// Array to hold FT(data)
	fftw_complex *FT_data_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Array to hold dFT(data)
	fftw_complex *deriv_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Compute Fourier transform of data
	ComputeFT(n, datainput, FT_data_c);
	
	
	// Multiply FT(data) by -ik to get derivatives in Fourier space
	// k = 2pi i / L
	double kdum;
	for(int i = 0; i < n; i++){
		kdum = 2.0 * PI * double(i) / ( h * n );
		// Compute first derivative in Fourier space
		// NOTE: Divide by "n" to normalise the Fourier coefficients
		deriv_FT_data[i][0] = -  kdum * FT_data_c[i][1] / n;
		deriv_FT_data[i][1] =    kdum * FT_data_c[i][0] / n;
		
	}
	
	// Compute   IFT( - ik FT( data ) )
	ComputeiFT(n, deriv_FT_data, datainput_derivFT);
	
	// Deallocate memory
	fftw_free(deriv_FT_data);
	fftw_free(FT_data_c);
	
	
} // END FFT_Deriv()

void FFT_DDeriv(int n, double h, double *datainput, double *datainput_dderivFT){

	// Function to return the second derivative
	// of real data, computed via FFT

	// Array to hold FT(data)
	fftw_complex *FT_data_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);	
	// Array to hold ddFT(data)
	fftw_complex *dderiv_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Compute Fourier transform of data
	ComputeFT(n, datainput, FT_data_c);
	
	// Multiply FT(data) -k^2 to get derivatives in Fourier space
	// k = 2pi i / L
	double kdum2;
	for(int i = 0; i < n; i++){
		kdum2 = pow( 2.0 * PI * double(i) / ( h * n ) , 2.0 );
		// Compute second derivative in Fourier space
		// NOTE:  Divide by "n" to normalise the Fourier coefficients
		for(int c = 0; c < 2; c++){
			dderiv_FT_data[i][c] = - kdum2 * FT_data_c[i][c] / n;
		}
	}
	
	// Compute   IFT( - k^2 FT( data ) )
	ComputeiFT(n, dderiv_FT_data, datainput_dderivFT);
	
	// Deallocate memory
	fftw_free(dderiv_FT_data);
	fftw_free(FT_data_c);
	
} // END FFT_Deriv()




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF