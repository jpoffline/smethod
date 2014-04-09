#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>


#include "inireader.h"

#include <boost/filesystem.hpp> // Used for creating the output filenames
#include <sstream> // Used for manipulating filenames


using namespace std;

const double pi = 4.0*atan(1.0);



void initialize_source_pot(int test_number, double h);
void iterate(int methodID);
void gauss_seidel(int step);
void successive_over_relaxation(int step);
double compute_error(int step);
void output_finalconfig(int methodID);
void checkdirexists(string dir);

