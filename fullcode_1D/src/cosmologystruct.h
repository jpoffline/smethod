
// cosmologystruct.h

// This contains all cosmology related stuff.
// It is held inside the field struct;
// and any cosmology related info can be obatined with the command
// field->cosmology.H0 (for today's Hubble value)
// field->cosmology.a (for a at the given eta)

#ifndef STRUCTCOSM_H
#define STRUCTCOSM_H

#include "fieldstruct.h"

struct COSM{
	
	// todays Hubble radius
	double H0;
	
	// length of the circle
	double L;	
	
	// \hbar
	double hbar; 
	
	// This is the value of superconformal time.
	// Computed as a function of the time-step
	double eta;
	
	// Scale factor: this variable is computed every time step,
	// as a function of eta
	double a;
	
	// Hubble: H(a)
	double H;
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set the background cosmology info,
	// at a given time-step.
	
	void SetBGcosmology(struct GRIDINFO *grid, struct FIELDCONTAINER *field){
		
		// Based on the time-step number set superconformal time, eta (time-step),
		field->cosmology.Seteta(grid, field);
		
		// then compute scale factor, a(eta),
		field->cosmology.Seta(field);
		
		// and finally Hubble, H(a).
		field->cosmology.SetH(field);
		
	} // END SetBGcosmology()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set the superconformal time, at a given time-step.
	
	void Seteta(struct GRIDINFO *grid, struct FIELDCONTAINER *field){
		
		// How do we set the "eta" (superconformal time) from
		// the time-step number and time-step size?
		// The eta should be chosen such that a=1, for eta_today
				
		// eta = (time-step - total # time-steps) * delta t
		// NOTE: the solver stop at time-step = ntimsteps - 1
		// and so the last value of eta_last = -ht.
		// If eta = 0 is reached, then a blows up to infinity.
		field->cosmology.eta = ( grid->t - grid->ntimsteps ) * grid->ht;
		
	} // END SetEta()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set the scale factor, a, at a given eta.
	
	void Seta( struct FIELDCONTAINER *field ){
		// a = 4 / ( eta^2 H0^2 )
		field->cosmology.a = 4.0 / ( pow( field->cosmology.eta  * field->cosmology.H0, 2.0 ) );
		
	} // END geta()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set Hubble, H(a)
	
	void SetH( struct FIELDCONTAINER *field ){
		// H^2 = H0^2 / a^3
		// H = sqrt(H^2)
		field->cosmology.H = field->cosmology.H0 * pow( field->cosmology.a, - 3.0 / 2.0 ) ;
		
	} // END getH()
	
};

#endif


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF