#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>

int main(void){

  int n = 2;
  int m = 4;
  int p = 2;


  int i, j, k;
  double *x1;
  fftw_complex *x2;

//////ALLOCATE MEM//////////////
  x1 = (double*) fftw_malloc(n*m*p*sizeof(double));
  x2 = (fftw_complex*) fftw_malloc(n*m*p*sizeof(fftw_complex));

/////SET UP FFT PLANS//////
   fftw_plan p1=fftw_plan_dft_r2c_3d(n,m,p, x1, x2,FFTW_MEASURE);
   fftw_plan p2=fftw_plan_dft_c2r_3d(n,m,p, x2, x1,FFTW_MEASURE);

//////CREATE INITIAL FUNCITON//////

printf("\n\nINITIAL 3D ARRAY\n");
for(k=0; k<p; k++){
  for (i=0; i < n; i++){
           for(j=0; j<m; j++){
                 x1[k+p*(j+m*i)]=sin(i+j*k*M_PI);
                printf("%lf\t", x1[k+p*(j+m*i)]);
        }printf("\n");
        }printf("\n");
        }

printf("\n\n");

//////FORWARD FFT//////////////////
fftw_execute(p1);
fftw_execute(p2);
printf("FORWARD FFTW RESULTS\n");
for(k=0; k<p; k++){
        for (i=0; i<n; i++){
                for(j=0; j<m; j++){
  printf("%g + %gi\t", creal(x2[k+p*(j*m*i)]),cimag(x2[k+p*(j+m*i)]));

}printf("\n");}printf("\n");
}
printf("\n\n");

//////BACKWARD FFT////////
//fftw_execute(p2);
printf("INVERSE FFTW RESULTS\n");
for(k=0; k<p; k++){
        for(i=0; i<n; i++){
                for(j=0; j<m; j++){
                printf("%lf\t", x1[k+p*(j*m*i)]/(n*p*m));
                }printf("\n");}printf("\n");
        }

//////DESTROY FFT PLANS///////
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);

//////FREE FFT ALLOCATED MEM///////
  fftw_free(x1);
  fftw_free(x2);
  return 0;
}