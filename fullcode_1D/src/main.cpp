
// KGsolve.cpp

// Main entry point into code.
// Setup to solve Klein-Gordon wave equation in 1D.
// User can specify the potential, equation of motion type,
// method to discretize derivatives.

#include "main.h"

// Can call executable with myParams.ini file
// - the default is params.ini

int main(int argc, char* argv[]){

// Start timing!
	boost::timer::cpu_timer myTimer;
	
	cout << endl;
	cout << "BEGIN" << endl;
	
// BEGIN: setup	

	// Two simple identifiers for writing parameters
	int ParamsPrintID_start = 0;
	int ParamsPrintID_end = 100;

	// The "field" struct contains everything to do with dynamical fields,
	// including the cosmological fields (like a & H)
	FIELDCONTAINER field;
	// The "params" struct contains all parameter choices etc,
	// and is mainly set by info read in from params.ini
	DATA params;
	// The "grid" struct contains all info about the grid the code
	// runs on (number of lattice sites, step sizes, etc)
	GRIDINFO grid;
	// The "poiss" struct contains all info about how to solve 
	// Poisson's equation
	POISS poiss;

	// Read in parameter files & populate "params" struct
	GetParams(argc, argv, &params);
	// Check the parameters for sanity
	CheckParams(&params);
	
	// If the parameters were ok, carry on:
	if( params.flag == 0){
	
		// Use params info to setup "grid",  "field", and "poiss" structs
		Setup(&params, &grid, &field, &poiss);
	
		// Print params to screen & logfile
		ofstream logout;
		logout.open(params.OutDir + params.RunID + "_log.dat");
		PrintParams(cout, &params, ParamsPrintID_start);	
		PrintParams(logout, &params, ParamsPrintID_start);	
		logout.close();
		
// END: setup	

// BEGIN: solving
		
		// Setup initial conditions.
		// This sets 
		// (a) cosmology
		// (b) field values
		InitialConditions(&params, &grid, &field);
		
		// Solve field equation.
		// Everything is done in here:
		// Run over time-steps,
		// run over space, solve Poisson's equation (if requested),
		// output to file,...
		SolveFieldEquation(&params, &grid, &field, &poiss);
		
		// Deallocate memory
		field.CleanField(&field);
		poiss.CleanPoiss(&poiss);
		
// END: solving

// BEGIN: feedback

// Stop timer
		myTimer.stop();
		params.TotalRunTime = myTimer.elapsed().wall / 1e6;
		
		// Finish off by writing elapsed time to the logfile.
		logout.open(params.OutDir + params.RunID + "_log.dat",std::ofstream::app);
		PrintParams(cout, &params, ParamsPrintID_end);	
		PrintParams(logout, &params, ParamsPrintID_end);	
		logout.close();
		
// END: feedback
		
	} // END if( params.flag == 0){}

	
}// end main()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF