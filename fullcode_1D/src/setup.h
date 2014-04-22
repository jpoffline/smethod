// setup.h

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


void SetupGrid(struct GRIDINFO *grid, struct DATA *params);
void SetupField(struct DATA *params, struct FIELDCONTAINER *field);
void SetupPoisson(struct DATA *params, struct POISS *poiss);
void checkdirexists(string dir);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF