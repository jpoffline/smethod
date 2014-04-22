
// poisssolve.h

#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>

using namespace std;

#include "gridstruct.h"
#include "datastruct.h"
#include "fieldstruct.h"

#include <fftw3.h>

void SolvePoisson(struct POISS *poiss);
void SolvePoisson_relax(struct POISS *poiss);
void SolvePoisson_FFT(struct POISS *poiss);

void ComputeFT(int n, double *datainput, fftw_complex *datainput_FT);
void ComputeiFT(int n, fftw_complex *datainput, double *datainput_iFT);



#define _USE_MATH_DEFINES

#define PI M_PI



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF