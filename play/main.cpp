


#include "main.h"

struct ARRAY{
    int ent1;
    int ent2;
};

int myarray_i(int i, int j, int k, int *max);
int doonstruct(ARRAY &arr);



int main(){
    
    int imax=10;
    int jmax=8;
    int kmax=3;
    int totlength=imax*jmax*kmax;
    double *myarray = new double[totlength];
    int *myarraym = new int[3];
    myarraym[0]=imax;
    myarraym[1]=jmax;
    myarraym[2]=kmax;
    
    for(int i=0;i<imax;i++){
        for(int j=0;j<jmax;j++){
            for(int k=0;k<kmax;k++){
                myarray[ myarray_i(i,j,k,myarraym)]=k;
            }
        }
    }
    cout << endl;
    cout << myarray[myarray_i(1,0,1,myarraym)] <<endl;
    cout << myarray[myarray_i(2,1,2,myarraym)] <<endl;
    cout << myarray[myarray_i(7,2,1,myarraym)] <<endl;
    cout << endl;
    delete myarray;
    delete myarraym;
    
    
    ARRAY arr1,arr2;
    arr1.ent1=4;
    arr1.ent2=8;
    arr2.ent1=5;
    
    int dum=doonstruct(arr1);
    cout <<dum << endl;
    
    cout << arr1.ent1 << " " << arr2.ent1 << endl;
    
    
}


int doonstruct(ARRAY &arr){
    return arr.ent1 * arr.ent2;
}


int myarray_i(int i, int j, int k, int *max){
    return i*max[1]*max[2]+j*max[2]+k;
}

