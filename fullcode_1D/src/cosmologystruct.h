
// cosmologystruct.h

// This contains all cosmology related stuff

#ifndef STRUCTCOSM_H
#define STRUCTCOSM_H


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
	
	void SetBGcosmology(struct GRIDINFO *grid, struct COSM *cosmology){
		
		// Based on the time-step number set superconformal time, eta (time-step),
		cosmology->Seteta(grid, cosmology);
		
		// then compute scale factor, a(eta),
		cosmology->Seta(cosmology);
		
		// and finally Hubble, H(a).
		cosmology->SetH(cosmology);
		
	} // END SetBGcosmology()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set the superconformal time, at a given time-step.
	
	void Seteta(struct GRIDINFO *grid, struct COSM *cosmology){
		
		// How do we set the "eta" (superconformal time) from
		// the time-step number and time-step size?
		// The eta should be chosen such that a=1, for eta_today
				
		// eta = (time-step - total # time-steps) * delta t
		// NOTE: the solver stop at time-step = ntimsteps - 1
		// and so the last value of eta_last = -ht.
		// If eta = 0 is reached, then a blows up to infinity.
		cosmology->eta = ( grid->t - grid->ntimsteps ) * grid->ht;
		
	} // END SetEta()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set the scale factor, a, at a given eta.
	
	void Seta( struct COSM *cosmology ){
		// a = 4 / ( eta^2 H0^2 )
		cosmology->a = 4.0 / ( pow( cosmology->eta  * cosmology->H0, 2.0 ) );
		
	} // END geta()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set Hubble, H(a)
	
	void SetH( struct COSM *cosmology ){
		// H^2 = H0^2 / a^3
		// H = sqrt(H^2)
		cosmology->H = cosmology->H0 * pow( cosmology->a, - 3.0 / 2.0 ) ;
		
	} // END getH()
	
};

#endif


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF