#include <cmath>
#include <float.h>
#include <fstream>
#include <string>

using namespace std;
#include <fftw3.h>
#define PI 3.1415926535897932
int main(){


	int numOfXGrid = 10;
	int numOfYGrid = 10;
	
	//These are the arrays in real space
	double **v;

	double *firstD_u;


	//arrays to store first order derivative
	double **v_x;


	

	//These are the arrays in fourier space
	fftw_complex **V;
	fftw_complex *temp_U;
	fftw_complex *firstD_U;


	//temporary files in storing derivatives from 1D to 2D
	fftw_complex **temp;

	//FFT transform


	fftw_plan plan_firstD;


	v = new double*[numOfXGrid];
	v_x = new double*[numOfXGrid];
	for(int i = 0; i < numOfXGrid; i++){
		v_x[i] = new double[numOfYGrid];
		v[i] = new double[numOfYGrid];
		for(int j = 0; j < numOfYGrid; j++){
			v_x[i][j] = 0;
			v[i][j]=sin(2*PI);

		}
	}
	
	V = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*numOfXGrid);
	temp_U = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*numOfXGrid*(numOfYGrid/2+1));
	firstD_U = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*numOfXGrid*(numOfYGrid/2+1));

	for(int i = 0; i < numOfXGrid; i++){
		V[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(numOfYGrid/2+1));
		for(int j = 0; j < numOfYGrid/2 + 1; j++){
			V[i][j][0] = 0;
			V[i][j][1] = 0;

		}
	}

	for(int i = 0; i < numOfXGrid*(numOfYGrid/2+1); i++){
		temp_U[i][0] = 0;
		temp_U[i][1] = 0;
		firstD_U[i][0] = 0;
		firstD_U[i][1] = 0;

	}

	temp = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*)*numOfXGrid);
	for(int i = 0; i < numOfXGrid; i++){
		temp[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(numOfYGrid/2+1));
	}
	/**=====================================================
	Initializing Plans
	======================================================*/
	

	



	/*=============================================
	         calculate dv/dx
	=============================================*/
	for(int i = 0; i < numOfYGrid/2 + 1; i++){
		for(int j = 0; j < numOfXGrid; j++){
			if(j < numOfXGrid/2){
				temp[j][i][0] = -2*PI*j/numOfXGrid*V[j][i][1];
				temp[j][i][1] = 2*PI*j/numOfXGrid*V[j][i][0];
			}
			if(j > numOfXGrid/2){
				temp[j][i][0] = -2*PI*(j-numOfXGrid)/numOfXGrid*V[j][i][1];
				temp[j][i][1] = 2*PI*(j-numOfXGrid)/numOfXGrid*V[j][i][0];
			}
		}
		temp[numOfXGrid/2][i][0] = 0;
		temp[numOfXGrid/2][i][1] = 0;
	}

	for(int i = 0; i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid/2+1; j++){
			firstD_U[i*(numOfYGrid/2+1)+j][0] = temp[i][j][0];
			firstD_U[i*(numOfYGrid/2+1)+j][1] = temp[i][j][1];
		}
	}
	plan_firstD = fftw_plan_dft_c2r_2d(numOfXGrid,numOfYGrid,firstD_U,firstD_u,FFTW_ESTIMATE);
	fftw_execute(plan_firstD);
	for(int i = 0;i < numOfXGrid; i++){
		for(int j = 0; j < numOfYGrid; j++){
			v_x[i][j] = firstD_u[i*numOfYGrid+j]/(numOfXGrid*numOfYGrid);
		}
	}
	

	
}// END main()