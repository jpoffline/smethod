
#include "readparams.h"

void GetParams(struct DATA *params){
    
    params->h = 0.1;
    params->ht = 0.01;
    params->imax = 10;
    params->jmax = 10;
    params->kmax = 10;
    params->cmax = 2;
    params->accuracy = 1E-4;
    params->derivsaccuracy = 2;
    
}