// Code to test accuracy of finite difference derivative schemes

#include <cmath>
#include <float.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/timer/timer.hpp>
using namespace std;
#include <fftw3.h>
#define _USE_MATH_DEFINES

#define PI M_PI
// Space step-size
double h = 0.1;
// Number of grid-points
int imax = 400;
// Function type (just for this testing purpose)
int fn_type = 3;
// How many waves in the box?
double wavn =  50.0;

double omega = wavn * 2.0 * PI / ( h * imax );
double xmax = h*imax;
double steep = 1.0;

double fn(double x){
	
	double ret;
	
	if(fn_type == 1) ret = sin( x * omega );
	if(fn_type == 2) ret = sin( x * omega ) + sin( x * 2.0 * omega );
	if(fn_type == 3) ret = tanh((0.25*xmax - x)/steep)+tanh((x-0.75*xmax)/steep);
	
	return ret;
	
}

double fn_x(double x){
	
	double ret;
	
	if(fn_type == 1) ret = omega*cos(x*omega);
	if(fn_type == 2) ret = omega*cos(x * omega) + 2*omega*cos(x*2.0*omega);

	return ret;

}

double fn_xx(double x){
	
	double ret;
	
	if(fn_type == 1) ret = -pow(omega,2.0)*sin(x*omega);
	if(fn_type == 2) ret = -pow(omega,2.0)*sin(x * omega) - pow(2*omega,2.0)*sin(x*2.0*omega);
	return ret;

}

void ComputeFT(int n, double *datainput, fftw_complex *datainput_FT){

	// This routine computes the Fourier transform of REAL data
	// Input: datainput[] -- dim n, data to be Fourier transformed
	// Output datainput_FT[] -- dim n, contains UNnormlised Fourier mode coefficients
	
	
	// Declare plan to compute FT(datainput)
	fftw_plan p; 

	// Construct plan
	p = fftw_plan_dft_r2c_1d(n, datainput, datainput_FT, FFTW_ESTIMATE);	
	
	// Execute plan
	fftw_execute(p); 
	
	// Delete plan
	fftw_destroy_plan(p);
	
	
} // END ComputeFT()

void ComputeiFT(int n, fftw_complex *datainput, double *datainput_iFT){

	// Declare plan to compute iFT(datainput)
	fftw_plan q;

	// Construct plan
	q = fftw_plan_dft_c2r_1d(n, datainput, datainput_iFT, FFTW_ESTIMATE);
	
	// Execute plan
	fftw_execute(q);
	
	// Delete plan
	fftw_destroy_plan(q);
	
} // END ComputeiFT()

void FFT_Deriv(int n, double *datainput, double *datainput_derivFT){

	// Function to return the first derivative
	// of real data, computed via FFT
	
	// Array to hold FT(data)
	fftw_complex *FT_data_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Array to hold dFT(data)
	fftw_complex *deriv_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Compute Fourier transform of data
	boost::timer::cpu_timer FFT_deriv_timer_1;
	ComputeFT(imax, datainput, FT_data_c);
	FFT_deriv_timer_1.stop();
	
	
	// Multiply FT(data) by -ik to get derivatives in Fourier space
	// k = 2pi i / L
	double kdum;
	for(int i = 0; i < n; i++){
		kdum = 2.0 * PI * double(i) / ( h * imax );
		// Compute first derivative in Fourier space
		// NOTE: Divide by "n" to normalise the Fourier coefficients
		deriv_FT_data[i][0] = -  kdum * FT_data_c[i][1] / n;
		deriv_FT_data[i][1] =    kdum * FT_data_c[i][0] / n;
		
	}
	
	// Compute   IFT( - ik FT( data ) )
	boost::timer::cpu_timer FFT_deriv_timer_2;
	ComputeiFT(imax, deriv_FT_data, datainput_derivFT);
	FFT_deriv_timer_2.stop();
	
	// Deallocate memory
	fftw_free(deriv_FT_data);
	fftw_free(FT_data_c);
	
	cout << endl;
	cout << "FT(data) :: " << FFT_deriv_timer_1.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "iFT(d_data) :: " << FFT_deriv_timer_2.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << endl;
	
} // END FFT_Deriv()

void FFT_DDeriv(int n, double *datainput, double *datainput_dderivFT){

	// Function to return the second derivative
	// of real data, computed via FFT

	// Array to hold FT(data)
	fftw_complex *FT_data_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);	
	// Array to hold ddFT(data)
	fftw_complex *dderiv_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	
	// Compute Fourier transform of data
	boost::timer::cpu_timer FFT_deriv_timer_31;
	ComputeFT(imax, datainput, FT_data_c);
	FFT_deriv_timer_31.stop();
	
	// Multiply FT(data) -k^2 to get derivatives in Fourier space
	// k = 2pi i / L
	double kdum2;
	for(int i = 0; i < n; i++){
		kdum2 = pow( 2.0 * PI * double(i) / ( h * imax ) , 2.0 );
		// Compute second derivative in Fourier space
		// NOTE:  Divide by "n" to normalise the Fourier coefficients
		for(int c = 0; c < 2; c++){
			dderiv_FT_data[i][c] = - kdum2 * FT_data_c[i][c] / n;
		}
	}
	
	// Compute   IFT( - k^2 FT( data ) )
	boost::timer::cpu_timer FFT_deriv_timer_32;
	ComputeiFT(imax, dderiv_FT_data, datainput_dderivFT);
	FFT_deriv_timer_32.stop();
	
	// Deallocate memory
	fftw_free(dderiv_FT_data);
	fftw_free(FT_data_c);
	cout << endl;
	cout << "FT(data) :: " << FFT_deriv_timer_31.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "iFT(dd_data) :: " << FFT_deriv_timer_32.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << endl;
	
} // END FFT_Deriv()

void FD_Deriv(int n, double *datainput, double *derivdata){
        
        // Routine to compute finite difference first derivative
        
        double  fip, fim;
        int     ip, im;
        for(int i = 0; i < n; i++){
			ip = i + 1;
			im = i - 1;
			if( ip == n ) ip = 0;
			if( im < 0 ) im = n - 1;

			fip = datainput[ip];
			fim = datainput[im];
			derivdata[i] = ( fip - fim ) / ( 2.0 * h );
        }
        
} // END FD_Deriv()

void FD_DDeriv(int n, double *datainput, double *derivdata){
        
        // Routine to compute finite difference second derivative        
        
        double  fip, fim, f0;
        int     ip, im;
        for(int i = 0; i < n; i++){
			ip = i + 1;
			im = i - 1;
			if( ip == n ) ip = 0;
			if( im < 0 ) im = n - 1;
			f0 = datainput[i];
			fip = datainput[ip];
			fim = datainput[im];
			derivdata[i] = ( fip + fim - 2*f0 ) / ( h * h );
        }
        
} // END FD_DDeriv()

int main(){
	
	cout << endl;
	cout << " -------------------------------------------" << endl;
	cout << "Code to test & compare accuracy of " << endl;
	cout << "finite difference and FFT derivative schemes" <<endl;
	cout << " > J Pearson, Durham, April 2014" << endl;
	cout << " -------------------------------------------" << endl;
	cout << "FD = finite difference" << endl;
	cout << "FFT = fast Fourier transform" << endl;
	cout << endl;
	string OutDir = "out/";
	
	double x;	
	double *data = new double[imax];
	double *fd_data = new double[imax];
	double *fdd_data = new double[imax];
	double *FFT_deriv = new double[imax];
	double *FFT_dderiv = new double[imax];	

	
	// Setup data
	for(int i = 0; i < imax; i++){
		x = i*h;
		data[i] = fn(x);
	}
	
	// Compute first derivative of data using finite difference
	boost::timer::cpu_timer fd_data_timer;
	FD_Deriv(imax, data, fd_data);
	fd_data_timer.stop();
	
	// Compute second derivative of data using finite difference
	boost::timer::cpu_timer fdd_data_timer;
	FD_DDeriv(imax, data, fdd_data);
	fdd_data_timer.stop();
	
	// Compute first derivative of data using FFT
	boost::timer::cpu_timer FFT_deriv_timer;
	FFT_Deriv(imax, data, FFT_deriv);
	FFT_deriv_timer.stop();
	
	// Compute second derivative of data using FFT
	boost::timer::cpu_timer FFT_dderiv_timer;
	FFT_DDeriv(imax, data, FFT_dderiv);
	FFT_dderiv_timer.stop();
	
	// Compute error (compared to exact expressions for derivatives)
	double error_FD_first=0.0;
	double error_FFT_first = 0.0;
	double error_FD_second=0.0;
	double error_FFT_second = 0.0;
		
	ofstream out;
	out.open(OutDir+"fd_test.dat");
	
	for(int i = 0; i < imax; i++){
	
		x = i * h;
		
		// Compute absolute error of each method at each location
		error_FD_first+= abs( fd_data[i] - fn_x(x) ) / imax;
		error_FFT_first+= abs( FFT_deriv[i] - fn_x(x) ) / imax;
		error_FD_second+= abs( fdd_data[i] - fn_xx(x) ) / imax;
		error_FFT_second+= abs( FFT_dderiv[i] - fn_xx(x) ) / imax;
		
		// Output coordinate & data
		out << x << " " << data[i] << " " ;
		// Output finite difference derivatives
		out << fd_data[i] << " " << fdd_data[i] << " " ;
		// Output FFT derivatives
		out << FFT_deriv[i] << " " << FFT_dderiv[i] << " " ;
		// Now output exact values: 
		out << fn(x) << " " << fn_x(x) << " " << fn_xx(x) << endl;		
		
	}

	out.close(); 
	cout << "Error on FD first derivative: " << error_FD_first << endl;
	cout << "Error on FD second derivative: " << error_FD_second << endl;
	cout << "Error on FFT first derivative: " << error_FFT_first << endl;
	cout << "Error on FFT second derivative: " << error_FFT_second << endl;
	cout << endl;
	
	cout << "FD first derivative compute time: " <<  fd_data_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "FD second derivative compute time: " <<  fdd_data_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "FFT first derivative compute time: " <<  FFT_deriv_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "FFT second derivative compute time: " <<  FFT_dderiv_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << endl;
	
	// Clean up
	delete data;
	delete fd_data;
	delete fdd_data;
	delete FFT_deriv;
	delete FFT_dderiv;
	
} // END main()