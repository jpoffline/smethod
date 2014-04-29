
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


void PrintParams(ostream& whereto, struct DATA *params, int ID);
void GetParams(struct DATA *params);
void SetupCosmology(struct DATA *params, struct COSM *cosmology);
void SetupGrid(struct GRIDINFO *grid, struct DATA *params);
void SetupField(struct DATA *params, struct FIELDCONTAINER *field);
void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct COSM *cosmology);
void SolveKG3D(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);
void GetDeriv(struct GRIDINFO *grid, struct FIELDCONTAINER *field);

void SetInitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct COSM *cosmology);





////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF