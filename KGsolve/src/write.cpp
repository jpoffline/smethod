
// write.cpp

#include "write.h"

void PrintParams(ostream& whereto,struct DATA *params, int ID){

	// Write this at the start of a run
	if( ID == 0 ){
		whereto << endl;
		whereto << "#------------------------------------------" << endl;
		whereto << "# Klein Gordon evolution engine" << endl;
		whereto << "# J. Pearson, Durham, April 2014" << endl;
		whereto << "#------------------------------------------" << endl;
		whereto << endl;
		whereto << "OutDir = " << params->OutDir << endl;
		whereto << "RunID = " << params->RunID << endl;
		whereto << endl;
		whereto << "h = " << params->h << endl;
		whereto << "ht = " <<  params->ht << endl;
		whereto << "imax = " << params->imax << endl;
		whereto << "jmax = " << params->jmax << endl;
		whereto << "kmax = " << params->kmax << endl;
		whereto << "cmax = " << params->cmax << endl;
		whereto << "ntimsteps = " << params->ntimsteps << endl;
		whereto << "filefreq = " << params->filefreq << endl;
		whereto << "accuracy = " << params->accuracy << endl;
		whereto << "derivsaccuracy = " << params->derivsaccuracy << endl;
		whereto << "inittype = " << params->inittype << endl;
		whereto << "evoltype = " << params->evoltype << endl;
		if( params->evoltype == 0){
			whereto << "# Gradient flow evolution" << endl;
		}
		if( params->evoltype == 1){
			whereto << "# 2nd order wave equation" << endl;
		}
		whereto << "pottype = " << params->pottype << endl;
		if( params->pottype == 0 ){
			whereto << "# Massive Scalar, mass = " << params->potparam1  << endl;
		}
		if( params->pottype == 1 ){
			whereto << "# Higgs potential" << endl;
		}
		whereto << endl;
    } // END if( ID == 0 ){}
    
    // Write this at the end of a run
    if( ID == 1 ){
	    whereto << endl;
    	if(params->flag == 0) whereto << "Completion without flag" << endl;
    	whereto << "Total run-time = " << params->TotalRunTime << " milliseconds" << endl;
    	whereto << endl;
    	whereto << "END" << endl;
    	whereto << endl;
    } // END if( ID == 1 ){}
    
} // END PrintParams()
