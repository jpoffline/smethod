
// kgsolve.h

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
#include "lapstencil.h"

void PrintParams(ostream& whereto, struct DATA *params);
void GetParams(struct DATA *params);
void SetupGrid(struct GRIDINFO *grid, struct DATA *params);
void SetupField(struct DATA *params, struct FIELDCONTAINER *field);
void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);
void SolveKG3D(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);
void GetDeriv(struct GRIDINFO *grid, struct FIELDCONTAINER *field);

