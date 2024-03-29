
// kgsolve.cpp

#include "kgsolve.h"

// This routine runs over time and space and solves the Klein-Gordon equation.
// Most of the calculations are done within the field struct

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


void SolveKG1D(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field ){
	
	// Create time-history struct:
	THIST timehistory;
	// Setup the time-history items, 
	timehistory.setup(&timehistory, params->thistfreq);	
	// and now zero the time-history counter.
	int th = 0;
	

	// This is the number used to identify field files output over time.
	// The value of "timefile" increments after each file is outputted.
	// Files are outputted with frequency params->filefreq
	int timefile = 100000;
	
	// This is a check to see if a file is being outputted at a given time-step.
	// fileout = 0 if no file being outputted;
	// fileout = 1 if a file is being outputted.
	int fileout = 0;
	
	// Some ofstreams used to dump data
	// ofstreams for field data:
	ofstream fieldout, finalfieldout;
	// ofstreams for numerical solutions to the Poisson equation:
	ofstream potout, finalpotout;
	
	// Open up time-history file
	timehistory.TimeHistoryFile.open( params->OutDir + params->RunID + "_timehistory.dat" );
	
	// Begin looping over time-steps
	for(int t = 0; t < params->ntimsteps; t++){
		
		// Set the time index for the field.
		// This function sets the given time-step number,
		// and the time-step identifiers (prev, now, next)
		// for the field values.
		grid->SetTime(t,grid);
		
		// Set the background cosmology values at this time-step.
		// This function sets 
		// (1) superconformal time, eta(time-step),
		// (2) scale factor, a(eta),
		// (3) Hubble, H(a).
		field->cosmology.SetBGcosmology(grid, field);

		// Dump stuff to file
		if( t % params->filefreq == 0 ){
		
			// (1) Change fileout flag 
			fileout = 1;
			
			// (2) Open up files
			fieldout.open( params->OutDir + params->RunID + "_f_" + to_string(timefile) + ".dat" );
			if( params->PoissSolnMethod != 0 ) potout.open( params->OutDir + params->RunID + "_V_" + to_string(timefile) + ".dat" );
			
			// (3) Increment filename flag (for next time)
			timefile++;
			
		}
		
		// Dump all info when t gets to the last time step.
		if( t == params->ntimsteps - 1 ){
		
			finalfieldout.open( params->OutDir + params->RunID + "_f_final.dat" );
			if( params->PoissSolnMethod != 0 ) finalpotout.open( params->OutDir + params->RunID + "_V_final.dat" );
			fileout = 10;
			
		}
				
		// Choose whether to compute Laplacian via FFT.
		// If we did, then use this in the equation of motion.
		// Slows code down by ~factor 3 or so, but more accurate than finite difference
		if( params->field_lap_type == 1 ) ComputeLaplacian_FFT(params, grid, field);
	
		// Solve Poisson's equation to get V.
		if( params->PoissSolnMethod != 0 ) {
			// Only solve it when timestep -modulo- PossSolveFreq = 0
			// since we dont expect Poissons equation to need solving
			// every time-step (although, it may need to be!).
			if( t % params->PossSolveFreq == 0 )
				SolvePoisson(params, grid, field);
		}
		
		// Run over the grid: compute EoM & update field
		//	- also, do any analysis (if required).
		for(int i = grid->imin; i < grid->imax; i++){
	
			// Set the location:
			// Current location i & forwards and backwards,
			// ip = i+1, im = i-1,
			// based on periodic boundary conditions.
			grid->GetPos(i,grid,0);
		
			// Do the bulk of solving the Klein-Gordon equation:
		
			// (1) Get the spatial derivatives of the field.
			// This either computes finite difference derivatives
			// or uses the FFT Laplacian computed just above this i-loop.
			// Which is used is chosen via "field_lap_type".
			field->GetDeriv(params, grid, field);
			
			// (2) Get derivative of the potential.
			// The user can specify which potential to use
			// via "pottype".
			field->equations.Getdpot(params, grid, field);
			
			// (3) Construct equation of motion.
			// Typically, this constructs E = \nabla^2\phi - V'(\phi),
			// but can also be used to construct e.g. the Schrodinger equation of motion;
			// which is chosen via "eomtype".
			field->equations.GetEoM(params, grid, field);
			
			// (4) Update value of the field.
			// This sets E = \dot{\phi} or E = \ddot{\phi}, depending on whether 
			// full evolution or gradient flow evolution is selected via "evoltype".
			// The field is then incremented
			field->UpdateField(params, grid, field);
			
			// Dump field values to file
			if( fileout > 0 ){				
				if( fileout == 1 ){
					field->WriteFieldData(fieldout,params,grid,field);
					if( params->PoissSolnMethod != 0 ) field->poiss.WritePoissData(i, potout, &field->poiss);
				}
				if( fileout == 10 ){
					field->WriteFieldData(finalfieldout, params, grid, field);							
					if( params->PoissSolnMethod != 0 ) field->poiss.WritePoissData(i, finalpotout, &field->poiss);
				}
			}
			
		} // END i-loop

		// Close the field-file (if it was open)
		if( fileout > 0 ){
			if( fileout == 1 ){
				fieldout.close();
				if( params->PoissSolnMethod != 0 )  potout.close();
				fileout = 0;
			}
			// Dump field configuration at the last time-step
			if( fileout == 10 ){
				finalfieldout.close();
				if( params->PoissSolnMethod != 0 )  finalpotout.close();
				fileout = 0;
			}
		}
		 
		// Construct time-history items.
		// These are set inside timehistorystruct.h
		timehistory.SetItems( th, grid, field, &timehistory );
		
		// Output info to screen
		if( t % params->screenfreq == 0 || t == params->ntimsteps - 1) {
			// Output time-step number & as % of everything to go
			cout << "(" << 100*t/(params->ntimsteps-1) << "%) ";
			// Dump time-history items to screen
			timehistory.TimeHistoryDump( cout, &timehistory, th );
		}
		
		// When the timestep number is the correct multiple of thistfreq,  
		// write the time-history file.
		if( t % params->thistfreq == 0 && th != 0){		
			timehistory.write( timehistory.TimeHistoryFile, &timehistory, th );
			th = -1;
		}
		// At the end of a run, make sure to dump the rest of the timehistory info
		if( t == params->ntimsteps - 1 && th != 0){		
			timehistory.write( timehistory.TimeHistoryFile, &timehistory, th );
			th = -1;
		}
		
		// increment counter for timehistory
		th++;
		
	} // END t-loop

	// Close the time-history output file,
	timehistory.TimeHistoryFile.close();
	// and deallocate the memory used by the time-history items.
	timehistory.CleanUp(&timehistory);
	
} // END SolveKG1D()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF