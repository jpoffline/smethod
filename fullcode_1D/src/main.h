
// main.h


#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>
#include "inireader.h"

using namespace std;

#include "gridstruct.h"
#include "datastruct.h"
#include "fieldstruct.h"
#include "poissstruct.h"
#include "cosmologystruct.h"

#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <boost/timer/timer.hpp>

// Functions called in main.cpp
void PrintParams(ostream& whereto, struct DATA *params, int ID);
void GetParams(int argc, char* argv[], struct DATA *params);
void CheckParams(struct DATA *params);
void Setup(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct POISS *poiss);
void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);
void SolveFieldEquation(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct POISS *poiss);


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF