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


void FFT_Deriv(int n, double h, double *datainput, double *datainput_derivFT);
void FFT_DDeriv(int n, double h, double *datainput, double *datainput_derivFT);
void ComputeFT(int n, double *datainput, fftw_complex *datainput_FT);
void ComputeiFT(int n, fftw_complex *datainput, double *datainput_iFT);


#define _USE_MATH_DEFINES

#define PI M_PI



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF