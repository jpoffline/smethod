
// initialconditions.h

#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>
#include <complex>

using namespace std;
typedef complex<double> dcmplx;

#define _USE_MATH_DEFINES

#define PI M_PI

#include "gridstruct.h"
#include "datastruct.h"
#include "fieldstruct.h" 
#include "cosmologystruct.h"



void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct COSM *cosmology);
void SetInitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct COSM *cosmology);





////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF