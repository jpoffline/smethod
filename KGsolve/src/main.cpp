
// KGsolve.cpp

// Main entry point into code.
// Setup to solve Klein-Gordon wave equation in 3D.
// User can specify the potential.

#include "main.h"

// Can call executable with myParams.ini file
// - the default is params.init

int main(int argc, char* argv[]){

// Start timing!
	boost::timer::cpu_timer myTimer;
	
	cout << endl;
	cout << "BEGIN" << endl;
	
// BEGIN: setup	

	FIELDCONTAINER field;
	DATA params;
	GRIDINFO grid;
	
	// Read in parameter files & populate "params"
	GetParams(argc,argv,&params);
	// Use info to setup "grid" and "field" struct
	SetupGrid(&grid, &params);
	SetupField(&params, &field);
	
	// Print params to screen & logfile
	ofstream logout;
	logout.open(params.OutDir + params.RunID + "_log.dat");
	PrintParams(cout, &params, 0);	
	PrintParams(logout, &params, 0);	
	logout.close();
		
// END: setup	

// BEGIN: solving

	// Setup initial conditions
	InitialConditions(&params, &grid, &field);
	// Solve field equation
	SolveKG3D(&params, &grid, &field);
	// Delete arrays
	field.CleanField(&field);
	
// END: solving

// BEGIN: feedback

	myTimer.stop();
	params.TotalRunTime = myTimer.elapsed().wall / 1e6;
	logout.open(params.OutDir + params.RunID + "_log.dat",std::ofstream::app);
	PrintParams(cout, &params, 1);	
	PrintParams(logout, &params, 1);	
	logout.close();
	
// END: feedback
	
}// end main()



