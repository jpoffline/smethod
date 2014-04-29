
// solve.cpp

// This is the driver file for solving the field equation.
// Called by main.cpp, and just calls the function to solve the Klein-Gordon equation, in kgsolv.cpp


#include "solve.h"

void SolveFieldEquation(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct POISS *poiss, struct COSM *cosmology){
	
	SolveKG1D(params, grid, field, poiss, cosmology);
	
} // END SolveFieldEquation()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF