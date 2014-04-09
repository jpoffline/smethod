#include "write.h"

void PrintParams(ostream& whereto,struct DATA *params){

    whereto << "h = " << params->h << endl;
    whereto << "ht = " <<  params->ht << endl;
    whereto << "imax = " << params->imax << endl;
    whereto << "jmax = " << params->jmax << endl;
    whereto << "cmax = " << params->cmax << endl;
    whereto << "accuracy = " << params->accuracy << endl;
    whereto << "derivs accuracy = " << params->derivsaccuracy << endl;
    
} // END PrintParams
