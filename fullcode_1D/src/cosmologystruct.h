// cosmologystruct.h

// This contains all cosmology related stuff

#ifndef STRUCTCOSM_H
#define STRUCTCOSM_H

struct COSM{
	
	// Newton's constant
	const double G = 1.0;
	// Initial density of matter
	const double rho0 = 1.0;
	// Length
	const double L = params->imax * params->h;	
	// \hbar
	const double hbar = 1.0; 
	// Hubble
	double H = 1.0;
	// Scale factor: this is computed every time step,
	// and its value can be called via geta().
	double a;
	
	// Function to set the scale factor,a, at a given time, \eta.
	void Seta(double eta, struct COSM *cosmology){
		// This is called at the start of a given time-step
		
		cosmology->a = 3.0 / ( pow( eta , 2.0 ) * 2.0 * PI * cosmology->G * cosmology->rho0 );
		
	} // END geta()
	
	// Function to return Hubble, H, at a given scale factor, a.
	double getH(double a, struct COSM *cosmology){
		
		return sqrt( 8 * PI * cosmology->G / 3.0 * cosmology->rho0 / pow( a, 3.0 ) );
		
	} // END getH()
	
};

#endif