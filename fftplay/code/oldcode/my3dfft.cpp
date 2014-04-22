#include <cmath>
#include <float.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/timer/timer.hpp>
using namespace std;
#include <fftw3.h>

#define PI 3.1415926535897932
#define spaceloop_1D_i for(i=0;i<imax; i++) 
#define spaceloop_1D_j for(j=0;j<jmax; j++) 
#define spaceloop_1D_k for(k=0;k<jmax; k++) 
#define spaceloop_3D for(i=0;i<imax; i++) for(j=0;j<jmax;j++) for(k=0;k<kmax;k++)

void FFT_MyData(int *dims, double *datainput, double* FT_of_data){

	int imax = dims[0];
	int jmax = dims[1];
	int kmax = dims[2];
	int maxlen = imax * jmax * kmax;
	int i, j, k, index;
	
	////////////////////////////////
	fftw_complex *FT_of_in;
	fftw_plan p;
	
	FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * maxlen);

	// Create a plan to compute Fourier transform
	p = fftw_plan_dft_r2c_3d(imax, jmax, kmax, datainput, FT_of_in, FFTW_ESTIMATE);	
	fftw_execute(p); 

	spaceloop_3D{
		index = (i * jmax + j ) * kmax + k ;
		FT_of_data[ index ] = FT_of_in[index][0]/maxlen;
		FT_of_data[  index + 1 ] = FT_of_in[index][1]/maxlen;
		if(abs(FT_of_in[index][0])>1E-5){
			cout <<"RE "<< i << " " << j << " " << k << " " << FT_of_in[index][0]<<endl;
		}
		if(abs(FT_of_in[index][1])>1E-5){
			cout <<"IM "<< i << " " << j << " " << k << " " << FT_of_in[index][1]<<endl;
		}
	}
	ofstream ou;
	ou.open("FT.dat");
	j=0;k=0;
	spaceloop_1D_i{
		index = (i * jmax + j ) * kmax + k ;
		ou << i << " " << FT_of_in[index][0] << " " << FT_of_in[index][1]<< endl;
	}
	ou.close();
	fftw_destroy_plan(p);
	fftw_free(FT_of_in);
		
} // END FFT_MyData()

void FFT3D_Deriv(int ID, int *dims, double *kx, double *ky, double *kz, double *datainput, double *derivdata){

	int imax = dims[0];
	int jmax = dims[1];
	int kmax = dims[2];
	int maxlen = imax * jmax * kmax;
	int i, j, k;
	
	////////////////////////////////
	fftw_complex *FT_of_in, *deriv_of_FT_of_in;
	fftw_plan p, q;
	
	FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * maxlen);
	deriv_of_FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * maxlen);
	// Create a plan to compute Fourier transform
	p = fftw_plan_dft_r2c_3d(imax, jmax, kmax, datainput, FT_of_in, FFTW_ESTIMATE);	
	fftw_execute(p); 
	fftw_destroy_plan(p);

	// Fourier space derivatives
	if(ID == 0){// Compute x-derivative
		spaceloop_3D{
			deriv_of_FT_of_in[(i*jmax+j)*kmax+k][0] =   kx[i] * FT_of_in[(i*jmax+j)*kmax+k][1] / maxlen;
			deriv_of_FT_of_in[(i*jmax+j)*kmax+k][1] = - kx[i] * FT_of_in[(i*jmax+j)*kmax+k][0] / maxlen;	
		}
	}
	if(ID == 1){ // Compute y-derivative
		spaceloop_3D{
			deriv_of_FT_of_in[(i*jmax+j)*kmax+k][0] =   ky[j] * FT_of_in[(i*jmax+j)*kmax+k][1] / maxlen;
			deriv_of_FT_of_in[(i*jmax+j)*kmax+k][1] = - ky[j] * FT_of_in[(i*jmax+j)*kmax+k][0] / maxlen;
		}
	}
	if(ID == 2){ // Compute z-derivative
		spaceloop_3D{
			deriv_of_FT_of_in[(i*jmax+j)*kmax+k][0] =   kz[k] * FT_of_in[(i*jmax+j)*kmax+k][1] / maxlen;
			deriv_of_FT_of_in[(i*jmax+j)*kmax+k][1] = - kz[k] * FT_of_in[(i*jmax+j)*kmax+k][0] / maxlen;
		}
	}
	
	// Create plan to compute inverse Fourier transform of the Fourier derivative
	q = fftw_plan_dft_c2r_3d(imax, jmax, kmax, deriv_of_FT_of_in, derivdata, FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);


	fftw_free(FT_of_in); 
	fftw_free(deriv_of_FT_of_in);
	////////////////////////////////
	

} // END FFT3D_Deriv()

void FFT3D_Laplacian(int *dims, double *kx, double *ky, double *kz, double *datainput, double *Lapdata){

	int imax = dims[0];
	int jmax = dims[1];
	int kmax = dims[2];
	int i, j, k, maxlen;
	maxlen = imax * jmax * kmax;
	double ksq;
	////////////////////////////////
	fftw_complex *FT_of_in, *lap_of_FT_of_in;
	fftw_plan p, q;
	
	FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * maxlen);
	lap_of_FT_of_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * maxlen);
	// Create a plan to compute Fourier transform
	p = fftw_plan_dft_r2c_3d(imax, jmax, kmax, datainput, FT_of_in, FFTW_ESTIMATE);	
	fftw_execute(p); 
	fftw_destroy_plan(p);

	// Fourier space Laplacian
	spaceloop_3D{
		ksq = pow(kx[i], 2.0) + pow(ky[j], 2.0) + pow(kz[k], 2.0);
		lap_of_FT_of_in[(i*jmax+j)*kmax+k][0] = - ksq * FT_of_in[(i*jmax+j)*kmax+k][0] / maxlen;
		lap_of_FT_of_in[(i*jmax+j)*kmax+k][1] = - ksq * FT_of_in[(i*jmax+j)*kmax+k][1] / maxlen;	
	}

	// Create plan to compute inverse Fourier transform on Fourier derivative
	q = fftw_plan_dft_c2r_3d(imax, jmax, kmax, lap_of_FT_of_in, Lapdata, FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);

	// Normalise derivatives
	//spaceloop_3D derivdata[(i*jmax+j)*kmax+k] *= 1.0/(imax*jmax*kmax);

	fftw_free(FT_of_in); 
	fftw_free(lap_of_FT_of_in);
	////////////////////////////////
	

} // END FFT3D_Deriv()

void FD_Laplacian(int *dims, double h, double *datainput, double *FD_lap_data){

	int imax = dims[0];
	int jmax = dims[1];
	int kmax = dims[2];
	double 	fip,fim,fjp,fjm,fkp,fkm,f0;
	int		i,j,k,ip,im,jp,jm,kp,km,index;
	double hx = h, hy = h, hz = h;

	spaceloop_3D{
	
		ip=i+1;
		im=i-1;
		jp=j+1;
		jm=j-1;
		kp=k+1;
		km=k-1;
		
		if( ip == imax ) ip = 0;
		if( im < 0 ) im = imax - 1;
		if( jp == jmax ) jp = 0;
		if( jm < 0 ) jm = jmax - 1;
		if( kp == kmax ) kp = 0;
		if( km < 0 ) km = kmax - 1;
		index = i * jmax * kmax + j * kmax + k;
		f0 = datainput[index];
		fip = datainput[ip * jmax * kmax + j * kmax + k];
		fim = datainput[im * jmax * kmax + j * kmax + k];
		fjp = datainput[i * jmax * kmax + jp * kmax + k];
		fjm = datainput[i * jmax * kmax + jm * kmax + k];
		fkp = datainput[i * jmax * kmax + j * kmax + kp];
		fkm = datainput[i * jmax * kmax + j * kmax + km];
		FD_lap_data[index] = (fip+fim-2.0*f0)/(hx*hx)
							+(fjp+fjm-2.0*f0)/(hy*hy)
							+(fkp+fkm-2.0*f0)/(hz*hz);

	}


} // END FD_Laplacian()

int main(){
	// Start timer
	boost::timer::cpu_timer myTimer;
	
	int imax = 100;
	int jmax = imax;
	int kmax = imax;
	double omx = 4.0;
	double omy = 0.0;
	double omz = 0.0;
	
	double h = 0.1;
	
	int ic_type = 1;
	
	
	int *dims = new int[3];
	dims[0] = imax;
	dims[1] = jmax;
	dims[2] = kmax;

	int i,j,k,index,maxlen;
	double x,y,z,val;
	
	cout << endl;
	cout << "3D FFT derivatives" << endl;
	cout << "Box-size: " << imax << "x"<< jmax << "x" << kmax << endl;
	cout << "omega_x = " << omx << endl;
	cout << "omega_y = " << omy << endl;
	cout << "omega_z = " << omz << endl;	
	cout << endl;

	maxlen = imax * jmax * kmax;
	double *kx = new double[imax];
	double *ky = new double[jmax];
	double *kz = new double[kmax];
	double *Data = (double*)fftw_malloc(imax*jmax*kmax*sizeof(double));
	double *Deriv_Data_x = (double*)fftw_malloc(maxlen*sizeof(double));
	double *Deriv_Data_y = (double*)fftw_malloc(maxlen*sizeof(double));
	double *Deriv_Data_z = (double*)fftw_malloc(maxlen*sizeof(double));	
	double *Deriv_Data_xx = (double*)fftw_malloc(maxlen*sizeof(double));
	double *Deriv_Data_yy = (double*)fftw_malloc(maxlen*sizeof(double));
	double *Deriv_Data_zz = (double*)fftw_malloc(maxlen*sizeof(double));	
	double *FD_lap_data = (double*)fftw_malloc(maxlen*sizeof(double));	
	double *FFT_Lapdata = (double*)fftw_malloc(maxlen*sizeof(double));	
	double *FT_of_data = (double*)fftw_malloc((2*maxlen+1)*sizeof(double));
	
	// setup wave-vectors
	/*
	spaceloop_1D_i kx[i] = double( i ) / ( h * imax );
	spaceloop_1D_j ky[j] = double( j ) / ( h * jmax );
	spaceloop_1D_k kz[k] = double( k ) / ( h * kmax );
	*/
	spaceloop_1D_i {
		kx[i] = (double(i)*h - 0.5*double(imax)*h);
		kx[i]=double(i)/(h*imax);
	}
	spaceloop_1D_j ky[j] = double( j ) / ( h * jmax );
	spaceloop_1D_k kz[k] = double( k ) / ( h * kmax );
	
	// Setup data to Fourier transform
	spaceloop_3D{
		index = i * jmax * kmax + j * kmax + k ;
		x = i*h; 
		y = j*h; 
		z = k*h;
		val = x*omx*2.0*PI/(h*imax) + y*omy*2.0*PI/(h*jmax) + z*omz*2.0*PI/(h*kmax);
		if( ic_type == 1 ){
			Data[ index ] = sin( val ) ;
		}
		if( ic_type == 2 ){
			Data[ index ] = sin( omx * kx[i] + omy * ky[j] + omz * kz[k] ) 
							+ sin( 2*omx * kx[i] + 3*omy * ky[j] + 2*omz * kz[k] );
		}
	}  	
	
	// Output slice of data
	ofstream outdum;
	outdum.open("ic_2d.dat");
	k = 0;
	spaceloop_1D_i{
		x = i*h;
		spaceloop_1D_j{
			index = i * jmax * kmax + j * kmax + k ;
			y = j*h;
			outdum << x << " " << y << " " << Data[ index ] << endl;
		}
		outdum <<endl;
	}
	outdum.close();
	
	// Get FT of data
	FFT_MyData(dims, Data, FT_of_data);
	outdum.open("FTdata.dat");
	j = 0; k = 0;
	spaceloop_1D_i{
		index = i * jmax * kmax + j * kmax + k;
		outdum << i << " " << FT_of_data[index] << " " << FT_of_data[index+1] << endl;
	}
	outdum.close();
	
	// Compute derivatives via FFT:
	boost::timer::cpu_timer FFT_timer;
	// (a) df/dx
	FFT3D_Deriv(0, dims, kx, ky, kz, Data, Deriv_Data_x);
	// (b) df/dy
	FFT3D_Deriv(1, dims, kx, ky, kz, Data, Deriv_Data_y);
	// (c) df/dz
	FFT3D_Deriv(2, dims, kx, ky, kz, Data, Deriv_Data_z);
	
	// (d) d2f/dx2
	FFT3D_Deriv(0, dims, kx, ky, kz, Deriv_Data_x, Deriv_Data_xx);
	// (e) d2f/dy2
	FFT3D_Deriv(1, dims, kx, ky, kz, Deriv_Data_y, Deriv_Data_yy);
	// (f) d2f/dz2
	FFT3D_Deriv(2, dims, kx, ky, kz, Deriv_Data_z, Deriv_Data_zz);
	FFT_timer.stop();
	
	boost::timer::cpu_timer FFT_lap_timer;
	FFT3D_Laplacian(dims, kx, ky, kz, Data, FFT_Lapdata);
	FFT_lap_timer.stop();
	
	boost::timer::cpu_timer FD_timer;
	// Compute Laplacian via finite difference
	FD_Laplacian(dims, h, Data, FD_lap_data);
	FD_timer.stop();
	
	// Output data & derivatives to file
	ofstream outfile;
	outfile.open("deriv_x.dat");
	j = 0; k = 0;
	spaceloop_1D_i{
		x = i * h;
		index = i * jmax * kmax + j * kmax + k;
		outfile << x << " " << Data[index] << " " << Deriv_Data_x[index] << " " << Deriv_Data_xx[index]  << endl;
	}
	outfile.close();
	
	outfile.open("deriv_y.dat");
	i = 0; k = 0;
	spaceloop_1D_j{
		y = j*h;
		index = i * jmax * kmax + j * kmax + k;
		outfile << y << " " << Data[index] << " " << Deriv_Data_y[index] << " "  << Deriv_Data_yy[index] << endl;
	}
	outfile.close();
	
	outfile.open("deriv_z.dat");
	i = 0; j = 0;
	spaceloop_1D_k{
		z = k*h;
		index = i * jmax * kmax + j * kmax + k;
		outfile << z << " " << Data[index] << " " << Deriv_Data_z[index]<< " " << Deriv_Data_zz[index] << endl;
	}
	outfile.close();
	
	outfile.open("deriv_lap.dat");
	j = 0; k = 0;
	spaceloop_1D_i{
		x = i*h;
		index = i * jmax * kmax + j * kmax + k;
		outfile << x << " " << Data[index] << " ";
		outfile << Deriv_Data_xx[index]+Deriv_Data_yy[index]+Deriv_Data_zz[index] << " " ;
		outfile << FFT_Lapdata[index] << " " ;
		outfile << FD_lap_data[index] << endl;
	}
	outfile.close();
	
	delete kx;
	delete ky;
	delete kz;
	delete dims;
	delete Deriv_Data_x;
	delete Deriv_Data_y;
	delete Deriv_Data_z;
	delete Deriv_Data_xx;
	delete Deriv_Data_yy;
	delete Deriv_Data_zz;
	delete FD_lap_data;
	delete FFT_Lapdata;
	delete FT_of_data;
	fftw_free(Data);
	
	/// Stop the timer & print elapsed time to screen
	myTimer.stop();
	double TotalRunTime = myTimer.elapsed().wall / 1e6;
	cout << "FFT Laplacian (method 1) compute time: " <<  FFT_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "FFT Laplacian (method 2) compute time: " <<  FFT_lap_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "FD Laplacian compute time: " <<  FD_timer.elapsed().wall / 1e6 << " milliseconds" << endl;
	cout << "total runtime = " << TotalRunTime << " milliseconds" << endl;
	cout << endl;
	
} // END main()

