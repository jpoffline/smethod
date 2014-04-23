
// solve.h

#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>

using namespace std;

#include "gridstruct.h"
#include "datastruct.h"
#include "fieldstruct.h"
#include "timehistorystruct.h"
#include "poissstruct.h"


void SolveKG1D(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct POISS *poiss);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF