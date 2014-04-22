#include <cmath>
#include <float.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
#include <fftw3.h>
#define PI 3.1415926535897932
#define spaceloop_1D for(int i=0; i<n; i++)

void FFT1D_Deriv(int n, double *datainput, double *derivdata, double *datainput_FT){

	// Setup wavenumbers
	double *kx = new double[n];
	spaceloop_1D{
		kx[i] = 2.0 * PI * i / n;
	}
	
	fftw_complex *FT_of_in, *deriv_of_FT_of_in;
	fftw_plan p, q;
	
	double *deriv_of_in = new double[n];
	FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	deriv_of_FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
		
	// Create a plan to compute Fourier transform
	p = fftw_plan_dft_r2c_1d(n, datainput, FT_of_in, FFTW_ESTIMATE);	
	fftw_execute(p); 
	fftw_destroy_plan(p);

	// Fourier space derivative
	spaceloop_1D{
		deriv_of_FT_of_in[i][0] = -  kx[i] * FT_of_in[i][1];
		deriv_of_FT_of_in[i][1] =    kx[i] * FT_of_in[i][0];
		for(int c = 0; c < 2; c++){
			datainput_FT[ 2 * i + c] = FT_of_in[i][c]/n;
		}
	}
	
	// Create plan to compute inverse Fourier transform
	q = fftw_plan_dft_c2r_1d(n, deriv_of_FT_of_in, derivdata, FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);

	delete deriv_of_in;
	delete kx;

	fftw_free(FT_of_in); 
	fftw_free(deriv_of_FT_of_in);


} // END FFT1D_Deriv()

void FFT1D_Laplacian(int n, double *datainput, double *LapData){

	// Setup wavenumbers
	double *kx = new double[n];
	spaceloop_1D{
		kx[i] = 2.0 * PI * i / n;
	}
	
	fftw_complex *FT_of_in, *lap_of_FT_of_in;
	fftw_plan p, q;
	
	double *deriv_of_in = new double[n];
	FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	lap_of_FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
		
	// Create a plan to compute Fourier transform
	p = fftw_plan_dft_r2c_1d(n, datainput, FT_of_in, FFTW_ESTIMATE);	
	fftw_execute(p); 
	fftw_destroy_plan(p);

	// Fourier space derivative
	spaceloop_1D{
		lap_of_FT_of_in[i][0] = - kx[i] * kx[i] * FT_of_in[i][0];
		lap_of_FT_of_in[i][1] = - kx[i] * kx[i] * FT_of_in[i][1];
	}
	
	// Create plan to compute inverse Fourier transform
	q = fftw_plan_dft_c2r_1d(n, lap_of_FT_of_in, LapData, FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);

	delete deriv_of_in;
	delete kx;

	fftw_free(FT_of_in); 
	fftw_free(lap_of_FT_of_in);


} // END FFT1D_Deriv()

// Finite difference derivative
void FD_Deriv(int n, double *datainput, double *derivdata){
	
	double 	fip, fim;
	int		ip, im;
	double h = 1.0/n;
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

double fn(double pos, double om, double totl){
	return sin( 2.0 * om *  PI * pos / totl ) ;
}
double dfn(double pos, double om){
	return 2.0 * om *  PI * cos( 2.0 * om *  PI * pos ) ;
}
double ddfn(double pos, double om){
	return - pow( 2.0 * om *  PI, 2.0) * sin( 2.0 * om *  PI * pos ) ;
}

int main(){

	string OutDir = "out/";

	int n = 400;
	double om = 5;
	double *datainput = new double[n];
	double *data_deriv_first = new double[n];
	double *data_deriv_fd = new double[n];
	double *data_deriv_second = new double[n];	
	double *data_FT = new double[2*n];
	double *deriv_FT = new double[2*n];
	double *LapData = new double[n];
	double h = 0.1;
	double TotPhysLength = double(n)*h;
	
	double pos;
	// Setup cos as input
	for(int i = 0; i < n; i ++){
		pos = double( i ) * h;
		datainput[i] = fn( pos, om , TotPhysLength );
	}
	
	// Compute finite difference derivative
	FD_Deriv(n, datainput, data_deriv_fd);
	// Compute FFT: first derivative
	FFT1D_Deriv(n, datainput, data_deriv_first, data_FT);
	// Compute FFT: second derivative
	FFT1D_Deriv(n, data_deriv_first, data_deriv_second, deriv_FT);

	// Compute Laplacian via FFT
	FFT1D_Laplacian(n, datainput, LapData);
	
	ofstream output;
	output.open(OutDir+"outdat.dat");
	double FFT_1std_error = 0.0;
	double FFT_2ndd_error = 0.0;
	double FFD_1d_error = 0.0;
	double exact_1, exact_2;
	
	spaceloop_1D{
		pos = double( i ) / double( n );
		exact_1 = dfn( pos, om );
		exact_2 = ddfn( pos, om );
		output << i << " " << datainput[i] << " " ;
		output << data_deriv_first[i]<< " " << data_deriv_fd[i] << " " ;
		output << data_deriv_second[i] << " " << LapData[i] << " " ;
		// Output analytic derivatives to compare
		output << exact_1 << " " ;
		output << exact_2 ;
		output << endl;
		FFT_1std_error += abs( (exact_1 - data_deriv_first[i]) );
		FFD_1d_error += abs( (exact_1 - data_deriv_fd[i]) );
		FFT_2ndd_error += abs( (exact_2 - data_deriv_second[i]) );
	}
	output.close();
	
	cout << "FFT 1st derivative error: " << FFT_1std_error << endl;
	cout << "Finite difference error: " << FFD_1d_error << endl;
	cout << "FFT 2nd derivative error: " << FFT_2ndd_error << endl;
		
	ofstream outputFT;
	outputFT.open(OutDir+"out_FourierModes.dat");
	double kx;
	
	spaceloop_1D{
		kx = 2.0 * PI * i / n;
		outputFT << kx << " " ;
		outputFT << data_FT[i] << " " << data_FT[i+1]  << " " ;
		outputFT << deriv_FT[i] << " " << deriv_FT[i+1]  << " " ;
		outputFT << endl;
	}
	outputFT.close();

	delete datainput;
	delete data_deriv_first;
	delete data_deriv_second;	
	delete data_deriv_fd;
	delete data_FT;
	delete deriv_FT;
	delete LapData;
	
} // END main()