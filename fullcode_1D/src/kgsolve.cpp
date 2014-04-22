
// kgsolve.cpp

#include "kgsolve.h"

// This routine runs over time and space and solves the Klein-Gordon equation.
// Most of the calculations are done within the field struct

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


void SolveKG3D(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){

	// Create time-history struct
	THIST timehistory;
	
	
	timehistory.SetFieldValDump(&timehistory);
	
	
	
	// This is the number used to identify field files output over time.
	// The value of "timefile" increments after each file is outputted.
	// Files are outputted with frequency params->filefreq
	int timefile = 100000;
	
	// This is a check to see if a file is being outputted at a given time-step.
	// = 0 if no file being outputted; = 1 if a file is being outputted.
	int fileout = 0;
	
	ofstream fieldout,fieldx, finalfieldout;
	
	
	// Open up time-history file
	timehistory.writeout.open( params->OutDir + params->RunID + "_timehistory.dat" );
	
	
	
	
	// Begin looping over time-steps
	for(int t = 0; t < params->ntimsteps; t++){
		
		// Set the time index for the field
		grid->SetTime(t,grid);
	
		// Output info to screen
		if( t % params->screenfreq == 0 ) cout << "(" << 100*t/params->ntimsteps << "%) Time-step number: "<< t << endl;
		
		// Dump stuff to file
		if( t % params->filefreq == 0 ){
		
			// (1) Change fileout flag 
			fileout = 1;
			// (2) Open up files
			fieldout.open( params->OutDir + params->RunID + "_" + to_string(timefile) + ".dat" );
			// (3) Increment filename flag (for next time)
			timefile++;
			
		}
		
		if( t == params->ntimsteps - 1 ){
		
			finalfieldout.open( params->OutDir + params->RunID + "_final.dat" );
			fileout = 10;
			
		}
				
		// Compute Laplacian via FFT
		// Will use this in the equation of motion
		// Slows code down by ~factor 3 or so, but more accurate than finite difference
		if( params->field_lap_type == 1 ) ComputeLaplacian_FFT(params, grid, field);
	
	
		// Solve Poisson's equation to get V
		/*
		SolvePoisson(poiss);
		*/
	
		// Run over the grid: compute EoM & update field
		//	- also, do any analysis (if required)
		for(int i = grid->imin; i < grid->imax; i++){
	
			grid->GetPos(i,grid,0);
		
			// (1) Get the spatial derivatives of the field
			field->GetDeriv(params, grid, field);
			// (2) Get derivative of the potential
			field->Getdpot(params, grid, field);
			// (3) Construct equation of motion
			field->GetEoM(params, field);
			// (4) Update value of the field
			field->UpdateField(params, grid, field);
			
			// Dump field values to file				
			if( fileout == 1 ) field->WriteFieldData(fieldout,params,grid,field);
			if( fileout == 10) field->WriteFieldData(finalfieldout, params, grid, field);							
			
		} // END i-loop

		// Close the field-file (if it was open)
		if( fileout == 1 ){
		
			fieldout.close();
			fileout = 0;

		}
		
		if(fileout == 10){
		
			finalfieldout.close();
			fileout = 0;
		
		}
		
		// Construct time-history items & write to file
		if( t%params->thistfreq == 0){
		
			timehistory.timestep = t;
			
			timehistory.time = t * params->ht;
			
			// Mainly for debugging & getting intuition,
			// dump the value of the field at a given location
			// at coordinates (valA_X), specified in timhistory_struct.h
			for(int v = 0; v < timehistory.nFieldVals_thist; v++){
				timehistory.val[v] = field->vals[field->ind(0,grid->now,timehistory.val_i[v],grid,field)];
			}

			timehistory.write(&timehistory);

		}
		
		
		
	} // END t-loop

	timehistory.writeout.close();
	timehistory.CleanUp(&timehistory);
	
	
	
} // END SolveKG3D()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF