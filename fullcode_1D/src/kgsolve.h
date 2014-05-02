
// kgsolve.h

#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include "gridstruct.h"
#include "datastruct.h"
#include "fieldstruct.h"
#include "poissstruct.h"
#include "timehistorystruct.h"



void ComputeLaplacian_FFT(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);
void SolvePoisson(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF