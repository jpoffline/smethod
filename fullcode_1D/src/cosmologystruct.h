// cosmologystruct.h

// This contains all cosmology related stuff

#ifndef STRUCTCOSM_H
#define STRUCTCOSM_H
#include <cmath>
#define _USE_MATH_DEFINES
#define PI M_PI
struct COSM{
	// todays Hubble radius
double H0;
	// length of the circle
double L;	
	// \hbar
double hbar; 
	// Scale factor: this is computed every time step,
	// and its value can be called via geta().
	double a;
	
	// Function to set the scale factor,a, at a given time, \eta.
	void Seta(struct DATA *params, struct GRIDINFO *grid, struct COSM *cosmology){
		// This is called at the start of a given time-step
		// the eta should be chosen such that a=1, for eta_today
		//
		
		cosmology->a = 4.0 / ( pow( (grid->now - params->ntimsteps)  * cosmology->H0, 2.0 ) );
		
	} // END geta()
	
	// Function to return Hubble, H, at a given scale factor, a.
	double getH(double a, struct COSM *cosmology){
		
		return H0 * sqrt( 1.0 / pow( a, 3.0 ) );
		
	} // END getH()
	
};

#endif