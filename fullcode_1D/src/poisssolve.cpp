
// poisssolve.cpp
 
#include "poisssolve.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Entry point into solving Poisson's equation.
// This chooses whether to use FFT or relaxation methods

void SolvePoisson(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field ){
	
	// Routine to choose how to solve Poisson equation
	// > nabla^2 V = S
	// - The potential "V" can be accesed in the main body
	// 	 of the codes
	
	// Setup the source, S, of the Poisson equation
	field->poiss.SetupSource( params, grid, field );
	
	// Solve via FFT
	if( field->poiss.method == 1 ){
		SolvePoisson_FFT( field );
	}
	
	// Solve via relaxation methods
	if( field->poiss.method == 2 ){
		SolvePoisson_relax( field );
	}
	
	// Compute error on the numerical solution
	// to the Poisson equation.
	 field->poiss.poisserr = Poisson_error( field );

} // END SolvePoisson()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to solve Poisson's equation via FFT

void SolvePoisson_FFT( struct FIELDCONTAINER *field ){
	
	// Solve Poisson's equation via FFT
	
	// nabla^2 V = S,
	
	// k = 2pi i / L, k2 = k^2
	double k, k2;
	double softening = 0.05;
	softening *= 2 * PI / ( field->poiss.h * field->poiss.imax );
	
	// Extract size, for easy reading of the code
	int n = field->poiss.imax;
	
	// Allocate memory
	// (1a) Shat = FT(S) -- Fourier transform of source, S
	fftw_complex *Shat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);	
	
	// (1b) Vhat = FT(V) -- Fourier transform of potential, V
	fftw_complex *Vhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Compute Fourier transform of source, S
	ComputeFT(n, field->poiss.S, Shat);

	// Divide FT(S) -k^2, to get Vhat = - Shat / k^2
	for(int i = 0; i < n; i++){
		k =  2.0 * PI * double(i) / ( field->poiss.h * n );
		//k+= softening;
		k2 = pow( k , 2.0 );
		
		for(int c = 0; c < 2; c++){
		
		// NOTE:  Divide by "n" to normalise the Fourier coefficients	
			Vhat[i][c] =  - Shat[i][c] / k2 / n;
			
			// Set Fourier coefficient to zero at k = 0
			// to ameliorate the pole.
			if( i == 0 ) Vhat[i][c] = 0.0;
			// This actually needs a lot more thought: how to deal with the zero mode!
			
		} // END c-loop
		
	} // END i-loop

	// Compute V = iFT(Vhat)
	ComputeiFT(n, Vhat, field->poiss.V);

	// Deallocate memory
	fftw_free(Shat);
	fftw_free(Vhat);
	

	
} // END SolvePoisson_FFT()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Entry point into solving Poisson's equation via relaxation methods.
// This chooses whether to use Gauss-Seidel or SOR.
// This loops over the given method until the error is less than desired accuracy

void SolvePoisson_relax( struct FIELDCONTAINER *field ){

	// Caller routine for solving Poissons equation
	// via relaxation algorithms
	// > Gauss Siedel & SOR.

	// Initialise the relaxation step counter
    int step = 0;
    // Holds the error on the "solution" to the Poisson equation
    double error;
    
    // Set initial guess on Poisson equation
    field->poiss.SetStep(step, field);    
	
    for(int i = 0; i < field->poiss.imax; i++){
    	field->poiss.rV[ field->poiss.rVi(field->poiss.now, i, field) ] = 0.0;
    }

    while(true){
    
        // Set the step -- leapfrogs inside each type
        // of relaxation routine
        field->poiss.SetStep(step, field);
        
        // Choose which relaxation method to use
        
        if( field->poiss.relaxmethod == 1) gauss_seidel( field );
            
        if( field->poiss.relaxmethod == 2) successive_over_relaxation( field );
        
        // Increment relaxation step
        step++;
        
        // Compute error on current "solution" to the Poisson equation
        error = relaxerror( field );
        
        // If error on "solution" is smaller than desired accuracy, stop
        if( error < field->poiss.accuracy )
            break;

    } // END loop for taking relaxation steps
    
    field->poiss.poisserr = error;
    
	// Dump result from relaxation into poiss->V[i]
	for(int i = 0; i < field->poiss.imax; i++){
		field->poiss.V[i] = field->poiss.rV[field->poiss.rVi(field->poiss.next, i, field)];
	}

} // END SolvePoisson_relax()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to compute error on "solution" to Poisson's
// equation.

double Poisson_error( struct FIELDCONTAINER *field ){

	int im, ip;
	double lap, source, error = 0.0;
	
	for(int i = 0; i < field->poiss.imax; i++){
	
		ip = i + 1;
		im = i - 1;
		if( ip == field->poiss.imax ) ip = 0;
		if( im < 0 ) im = field->poiss.imax - 1;
		lap = ( field->poiss.V[ip] + field->poiss.V[im] - 2.0 * field->poiss.V[i] ) / field->poiss.h2;
		source = field->poiss.S[i];
		error += abs(lap - source);
		
	} // END i-loop
	
	error = error / field->poiss.imax;
	return error;

} // END Poisson_error()

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF