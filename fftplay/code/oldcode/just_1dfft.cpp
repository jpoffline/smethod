#include <cmath>
#include <float.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
#include <fftw3.h>
#define PI 3.1415926535897932


void ComputeFT(int n, double *inputdata, fftw_complex *FT_of_data){

	
	fftw_plan p;	

	p = fftw_plan_dft_r2c_1d(n, inputdata, FT_of_data, FFTW_ESTIMATE);	
	
	fftw_execute(p); 
	fftw_destroy_plan(p);
	ofstream of;
	of.open("of.dat");
	for(int i=0; i<n; i++){
		of << i << " " << FT_of_data[i][0] << " " << FT_of_data[i][1] << endl;
	}
	of.close();
}

int main(){


	double Nwaves = 2.0;
	int imax = 200;
	double *datainput = new double[imax];
	double x;
	double DomainLength = 1.0;
	double h = DomainLength/double(imax) ;
	fftw_complex *FT_of_data  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * imax);;
	
	ofstream ic;
	ic.open("ic.dat");
	for(int i = 0; i < imax; i ++){
		x = i*h;
		datainput[i] = sin(Nwaves * 2.0*PI*x/DomainLength);
		ic << x << " " << datainput[i] << endl;
	}
	ic.close();
	
	ComputeFT(imax, datainput, FT_of_data);

	

	fftw_free(FT_of_data);
	delete datainput;

} // END main()