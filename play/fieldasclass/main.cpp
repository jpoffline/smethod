#include <iostream>
#include <cmath>
#include <float.h>
#include <fstream>
#include <string>

using namespace std;

void DoWithMyStruct(struct DATA *dat);

struct DATA{
    
    double h;
    int imax;
    int jmax;
    int kmax;
    
    
};


int main(){

    DATA params;
   
    DoWithMyStruct( &params);
}

void DoWithMyStruct(struct DATA *params){
    params->h=0.1;
    params->imax=10;
    params->jmax=10;
    cout << params->h <<endl;
    
    cout << params->h << endl;
    cout << params->imax << endl;
    cout << params->jmax << endl;


}