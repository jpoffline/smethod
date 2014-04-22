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
double hx = 0.2;
double hy = hx;
double hz = hx;
// Number of grid-points


void Compute3DFT(int *dims, double *datainput, fftw_complex *datainput_FT){

	// This routine computes the Fourier transform of REAL data
	// Input: datainput[] -- dim n, data to be Fourier transformed
	// Output datainput_FT[] -- dim n, contains UNnormlised Fourier mode coefficients
	
	fftw_plan p;

	// Create a plan to compute Fourier transform
	p = fftw_plan_dft_r2c_3d(dims[0],dims[1], dims[2], datainput, datainput_FT, FFTW_ESTIMATE);			

	// Execute plan
	fftw_execute(p); 
	
	// Delete plan
	fftw_destroy_plan(p);
	
	
} // END Compute3DFT()

void Compute3DiFT(int *dims, fftw_complex *datainput, double *datainput_iFT){

	// Declare plan to compute iFT(datainput)
	fftw_plan q;

	// Construct plan
	q = fftw_plan_dft_c2r_3d(dims[0],dims[1], dims[2], datainput, datainput_iFT, FFTW_ESTIMATE);	
	
	// Execute plan
	fftw_execute(q);
	
	// Delete plan
	fftw_destroy_plan(q);
	
} // END Compute3DiFT()

void FFT_Deriv(int *dims, double *datainput, double *datainput_derivFT){

	// Function to return the first derivative
	// of real data, computed via FFT

	double kx,ky,kz,kx2,ky2,kz2;
	int imax = dims[0];
	int jmax = dims[1];
	int kmax = int(dims[2]/2)+1;
	int max = imax * jmax * kmax, index;

	// Array to hold FT(data)
	fftw_complex *FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);
	
	// Array to hold dFT(data)
	fftw_complex *deriv_x_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);
	fftw_complex *deriv_y_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);
	fftw_complex *deriv_z_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);
	fftw_complex *deriv_xx_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);
	fftw_complex *deriv_yy_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);			
	fftw_complex *deriv_zz_FT_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * max);				
	
	// Compute Fourier transform of data
	boost::timer::cpu_timer FFT_deriv_timer_1;
	Compute3DFT(dims, datainput, FT_data);
	FFT_deriv_timer_1.stop();
	
	
	// Multiply FT(data) by -ik to get derivatives in Fourier space
	// k = 2pi i / L
	 
	for(int i = 0; i < imax; i++){
		kx = 2.0 * PI * double(i) / ( hx * imax );
		kx2 = pow(kx, 2.0);
		for(int j = 0; j < jmax; j++){
			ky = 2.0 * PI * double(j) / ( hy * jmax );
			ky2 = pow(ky, 2.0);
			for(int k = 0; k < kmax; k++){
				kz = 2.0 * PI * double(k) / ( hz * kmax );
				kz2 = pow(kz, 2.0);
				index = (i * jmax + j) * kmax + k;
				deriv_x_FT_data[index][0] = - kx * FT_data[index][1] / max;
				deriv_x_FT_data[index][1] =   kx * FT_data[index][0] / max;
				
				deriv_y_FT_data[index][0] = - ky * FT_data[index][1] / max;
				deriv_y_FT_data[index][1] =   ky * FT_data[index][0] / max;
				
				deriv_z_FT_data[index][0] = - kz * FT_data[index][1] / max;
				deriv_z_FT_data[index][1] =   kz * FT_data[index][0] / max;
				
				for(int c = 0; c < 2; c++){
					deriv_xx_FT_data[index][c] = -  kx2 * FT_data[index][c] / max;
					deriv_yy_FT_data[index][c] = -  ky2 * FT_data[index][c] / max;				
					deriv_zz_FT_data[index][c] = -  kz2 * FT_data[index][c] / max;
				}
				
			}
		}
	}
	
	
	ofstream out;
	out.open("DUM.dat");
	int jj=0, kk=0;
	for(int i=0; i<imax; i++){
		index = (i * jmax + jj) * kmax + kk;
		out << i << " "  << FT_data[index][0]/max << " " << FT_data[index][1]/max << endl;
	}
	out.close();
	
	// Compute   IFT( - ik FT( data ) )
	boost::timer::cpu_timer FFT_deriv_timer_2;
	Compute3DiFT(dims, deriv_x_FT_data, datainput_derivFT);
	FFT_deriv_timer_2.stop();
	
	// Deallocate memory
	fftw_free(deriv_x_FT_data);
	fftw_free(deriv_y_FT_data);
	fftw_free(deriv_z_FT_data);
	fftw_free(deriv_xx_FT_data);
	fftw_free(deriv_yy_FT_data);
	fftw_free(deriv_zz_FT_data);					
	fftw_free(FT_data);
	
	cout << endl;
	cout << "FT(data) :: " << FFT_deriv_timer_1.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "iFT(d_data) :: " << FFT_deriv_timer_2.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << endl;
	
} // END FFT_Deriv()

void FD_derivs(int *dims, double *datainput, double *FD_derivs){

	// FD_derivs[nderivs*imax*jmax*kmax]

	int imax = dims[0];
	int jmax = dims[1];
	int kmax = dims[2];
	int nderivs = dims[3];
	double fip,fim,fjp,fjm,fkp,fkm,f0;
	int	 ip,im,jp,jm,kp,km,index;
	
	for(int i = 0; i < imax; i++){
		ip=i+1;
		im=i-1;
		if( ip == imax ) ip = 0;
		if( im < 0 ) im = imax - 1;
		for(int j=0; j < jmax; j++){
			jp=j+1;
			jm=j-1;
			if( jp == jmax ) jp = 0;
			if( jm < 0 ) jm = jmax - 1;
			for(int k=0; k < kmax; k++){
			
				kp=k+1;
				km=k-1;
				if( kp == kmax ) kp = 0;
				if( km < 0 ) km = kmax - 1;

				index = (i * jmax + j) * kmax + k;
				f0 = datainput[index];
				fip = datainput[(ip * jmax + j ) * kmax + k];
				fim = datainput[(im * jmax + j ) * kmax + k];
				fjp = datainput[(i * jmax + jp ) * kmax + k];
				fjm = datainput[(i * jmax + jm ) * kmax + k];
				fkp = datainput[(i * jmax + j ) * kmax + kp];
				fkm = datainput[(i * jmax + j ) * kmax + km];
				FD_derivs[nderivs*index+0] = (fip-fim)/(2.0*hx);
				FD_derivs[nderivs*index+1] = (fjp-fjm)/(2.0*hy);
				FD_derivs[nderivs*index+2] = (fkp-fkm)/(2.0*hz);
				FD_derivs[nderivs*index+3] = (fip+fim-2.0*f0)/(hx*hx);
				FD_derivs[nderivs*index+4] = (fjp+fjm-2.0*f0)/(hy*hy);
				FD_derivs[nderivs*index+5] = (fkp+fkm-2.0*f0)/(hz*hz);			
									
			} // END k-loop						
		} // END j-loop
	} // END i-loop


} // END FD_derivs()


int main(){
	
	cout << endl;
	cout << " -------------------------------------------" << endl;
	cout << "Code to test & compare accuracy of " << endl;
	cout << "3D finite difference and FFT derivative schemes" <<endl;
	cout << " > J Pearson, Durham, April 2014" << endl;
	cout << " -------------------------------------------" << endl;
	cout << "FD = finite difference" << endl;
	cout << "FFT = fast Fourier transform" << endl;
	cout << endl;
	string OutDir = "out/";
	
	int imax = 50;
	int jmax = 2;
	int kmax = 2;

	int nderivs = 6;
	int ndims = 4;
	double x, y, z;	
	int index, max = imax * jmax * kmax;
	
	int *dims = new int[ndims];
	double *data = new double[max];
	double *derivs_FD = new double[nderivs*max];
	double *derivs_FFT = new double[max];
	
	dims[0] = imax;
	dims[1] = jmax;
	dims[2] = kmax;
	dims[3] = nderivs;
	
	// Setup data
	double wavn_x = 2.0;
	double wavn_y = 0.0;
	double wavn_z = 0.0;		
	double omega_x = wavn_x * 2.0 * PI / ( hx * imax );
	double omega_y = wavn_y * 2.0 * PI / ( hy * jmax );
	double omega_z = wavn_z * 2.0 * PI / ( hz * kmax );
	
	for(int i = 0; i < imax; i++){
		x = i * hx;
		for(int j = 0; j < jmax; j++){
			y = j * hy;
			for(int k = 0; k < kmax; k++){
				z = k * hz;
				index = (i * jmax + j) * kmax + k;
				data[index] = sin( x * omega_x + y * omega_y + z * omega_z);
				
			} // END k-loop
		} // END j-loop
	} // END i-loop
	
	// Compute first derivative of data using finite difference
	boost::timer::cpu_timer FD_derivs_timer;
	FD_derivs(dims, data, derivs_FD);
	FD_derivs_timer.stop();
	
	// Compute first derivative of data using FFT
	boost::timer::cpu_timer FFT_derivs_timer;
	FFT_Deriv(dims, data, derivs_FFT);
	FFT_derivs_timer.stop();
	 
	// Dump field & derivatives into file
	// (along x-axis)	
	ofstream out;
	out.open(OutDir+"fd_3d_test.dat");
	int jj = 0, kk = 10;
	for(int i = 0; i < imax; i++){
		x = i * hx;
		index = (i * jmax + jj) * kmax + kk;
		out << x << " " << data[index] << " " ;
		for(int d = 0; d < nderivs; d++){
			out << derivs_FD[nderivs*index+d] << " " ;
		}
		out << derivs_FFT[index] << " " ;
		out <<endl;
	}
	out.close();
	
 	cout << endl;
 	cout << "Finite difference derivatives :: compute time = " << FD_derivs_timer.elapsed().wall/1E6 << " milliseconds" << endl;
	cout << "FFT derivatives :: compute time = " << FFT_derivs_timer.elapsed().wall/1E6 << " milliseconds" << endl;
 	cout << endl; 
 	
	// Clean up
	delete dims;
	delete data;
	delete derivs_FD;
	delete derivs_FFT;
	
} // END main()