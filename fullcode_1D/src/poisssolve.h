
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
#include "poissstruct.h"

#include <fftw3.h>

void SolvePoisson( struct FIELDCONTAINER *field );
void SolvePoisson_relax( struct FIELDCONTAINER *field );
void SolvePoisson_FFT( struct FIELDCONTAINER *field );

void ComputeFT(int n, double *datainput, fftw_complex *datainput_FT);
void ComputeiFT(int n, fftw_complex *datainput, double *datainput_iFT);
void gauss_seidel( struct FIELDCONTAINER *field );
void successive_over_relaxation( struct FIELDCONTAINER *field );
double relaxerror( struct FIELDCONTAINER *field );
double Poisson_error( struct FIELDCONTAINER *field );
#define _USE_MATH_DEFINES

#define PI M_PI



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF