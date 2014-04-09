#include "main.h"

// main.cpp


struct FIELDARRAY{
  
    int ncom;
    int ni;
    int nj;
    int nk;
    int I(int i, int j, int k, struct FIELDARRAY *F){
        return ( i * F->nj + j ) * F->nk + k;
    }
    
    double lap;
    double dpot;
    double eom;
    
    double *arr;
 
    double Update(struct FIELDARRAY *F, struct DATA *params){
        return F->eom * params->ht;
    }
 
    
} ;

struct

void messwiths(double *F);

int main(){


	// Struct holding all parameters
    DATA params;
    // Struct holding all field quantities
    FIELDARRAY field;
    GetParams(&params);
    
    PrintParams(cout, &params);
    
    
    field.ni = params.imax;
    field.nj = params.jmax;
    field.nk = params.kmax;
    field.arr = new double[field.ni*field.nj*field.nk];
    
    

    
    for(int i=0;i<field.ni;i++){
        for(int j=0;j<field.nj;j++){
            for(int k=0;k<field.nk;k++){
                
                field.lap=1.0;
                field.dpot=1.0;
                field.eom=field.lap-field.dpot;
                field.Update(&field,&params);
                
            }
        }
    }
    
    
    
    
    cout << "q"  << endl;
    cout << "q"  << endl;
    cout << "q"  << endl;

    delete field.arr;

   /*
    writetoscreen();
    
    doinitialconditions();
    
    evolve();
    
    writefinalconfig();
    */
    
} // end main()
