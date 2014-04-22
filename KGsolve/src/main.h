

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
#include "lapstencil.h"

#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <boost/timer/timer.hpp>

void PrintParams(ostream& whereto, struct DATA *params, struct LAPLACIANSTENCIL *stencil, int ID);
void GetParams(int argc, char* argv[], struct DATA *params);
void SetupGrid(struct GRIDINFO *grid, struct DATA *params);
void SetupField(struct DATA *params, struct FIELDCONTAINER *field);
void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field);
void SolveKG3D(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct LAPLACIANSTENCIL *stencil);
void GetDeriv(struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct LAPLACIANSTENCIL *stencil);
void CheckParams(struct DATA *params);
void SetupLaplacianStencil(struct DATA *params, struct LAPLACIANSTENCIL *stencil);
void WriteLog(struct DATA *params, double TotalRunTime);