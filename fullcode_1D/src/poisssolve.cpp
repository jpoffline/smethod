
// poisssolve.cpp
 
#include "poisssolve.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Entry point into solving Poisson's equation.
// This chooses whether to use FFT or relaxation methods

void SolvePoisson(struct POISS *poiss){
	
	// Routine to choose how to solve Poisson equation
	// > nabla^2 V = S
	// - The potential "V" can be accesed in the main body
	// 	 of the codes
	
	// Setup the source, S, of the Poisson equation
	poiss->SetupSource(params, grid, field, poiss);
	
	// Solve via FFT
	if(poiss->method == 1){
		SolvePoisson_FFT( poiss );
	}
	
	// Solve via relaxation methods
	if(poiss->method == 2){
		SolvePoisson_relax( poiss );
	}

} // END SolvePoisson()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to solve Poisson's equation via FFT

void SolvePoisson_FFT(struct POISS *poiss){
	
	// Solve Poisson's equation via FFT
	
	// nabla^2 V = S,
	
	// k = 2pi i / L, k2 = k^2
	double k2;
	
	// Extract size, for easy reading of the code
	int n = poiss->imax;
	
	// Allocate memory
	// (1a) Shat = FT(S) -- Fourier transform of source, S
	fftw_complex *Shat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);	
	
	// (1b) Vhat = FT(V) -- Fourier transform of potential, V
	fftw_complex *Vhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Compute Fourier transform of source, S
	ComputeFT(n, poiss->S, Shat);

	// Divide FT(S) -k^2, to get Vhat = - Shat / k^2
	for(int i = 0; i < n; i++){
	
		k2 = pow( 2.0 * PI * double(i) / ( poiss->h * n ) , 2.0 );
	
		for(int c = 0; c < 2; c++){
		
		// NOTE:  Divide by "n" to normalise the Fourier coefficients	
			Vhat[i][c] = - Shat[i][c] / kdum2 / n;
		
		} // END c-loop
		
	} // END i-loop

	// Compute V = iFT(Vhat)
	ComputeiFT(n, Vhat, poiss->V);

	// Deallocate memory
	fftw_free(Shat);
	fftw_free(Vhat);
	
	
} // END SolvePoisson_FFT()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Entry point into solving Poisson's equation via relaxation methods.
// This chooses whether to use Gauss-Seidel or SOR.
// This loops over the given method until the error is less than desired accuracy

void SolvePoisson_relax(struct POISS *poiss){

	// Caller routine for solving Poissons equation
	// via relaxation algorithms
	// > Gauss Siedel & SOR.

	// Initialise the relaxation step counter
    int step = 0;
    // Holds the error on the "solution" to the Poisson equation
    double error;
    
    while(true){
    
        // Set the step -- leapfrogs inside each type
        // of relaxation routine
        poiss->SetStep(step);
        
        // Choose which relaxation method to use
        if( poiss->relaxmethod == 1) gauss_seidel( poiss );
            
        if( poiss->relaxmethod == 2) successive_over_relaxation( poiss );
        
        // Increment relaxation step
        step++;
        
        // Compute error on current "solution" to the Poisson equation
        error = relaxerror( poiss );
        
        // If error on "solution" is smaller than desired accuracy, stop
        if( error < poiss->accuracy )
            break;
        
    } // END loop for taking relaxation steps
    
	// NEED TO DUMP RESULT INTO poiss->V[i]

} // END SolvePoisson_relax()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF