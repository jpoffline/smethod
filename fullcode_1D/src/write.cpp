
// write.cpp

// This contains all the info which is written to screen and to logfile

#include "write.h"

void PrintParams(ostream& whereto, struct DATA *params, int ID){

	// Write this at the start of a run
	if( ID == 0 ){
	
		whereto << endl;
		whereto << "#------------------------------------------" << endl;
		whereto << "# 1D Klein Gordon evolution engine" << endl;
		whereto << "# J. Pearson, Durham, April 2014" << endl;
		whereto << "#------------------------------------------" << endl;
		whereto << endl;
		whereto << "OutDir = " << params->OutDir << endl;
		whereto << "RunID = " << params->RunID << endl;
		whereto << endl;
		whereto << "h = " << params->h << endl;
		whereto << "ht = " <<  params->ht << endl;
		whereto << "imax = " << params->imax << endl;
		whereto << "cmax = " << params->cmax << endl;
		whereto << "ntimsteps = " << params->ntimsteps << endl;
		whereto << "filefreq = " << params->filefreq << endl;
		whereto << "inittype = " << params->inittype << endl;
		
		whereto << "H0 = " << params->H0 << endl;
		whereto << "hbar = " << params->hbar << endl;		
		
		whereto << "eomtype = " << params->eomtype << endl;
		if( params->eomtype == 0 ){
			whereto << " # Wave equation type EoM: nabla^2phi - dV/dphi" << endl;
		}
		if( params->eomtype == 1 ){
			whereto << " # Schrodinger equation type EoM" << endl;
		}
		
		
		whereto << "evoltype = " << params->evoltype << endl;
		if( params->evoltype == 0){
			whereto << " # Gradient flow evolution" << endl;
		}
		if( params->evoltype == 1){
			whereto << " # 2nd order wave equation" << endl;
		}
		
		whereto << "field_lap_type = " << params->field_lap_type << endl;
		if(params->field_lap_type == 0){
			whereto << " # Finite difference Laplacian" << endl;
			whereto << "derivsaccuracy = " << params->derivsaccuracy << endl;
		}
		if(params->field_lap_type == 1){
			whereto << " # FFT Laplacian" << endl;
		}
		
		whereto << "pottype = " << params->pottype << endl;
		if( params->pottype == 0 ){
			whereto << " # Massive Scalar, mass = " << params->potparam1  << endl;
		}
		if( params->pottype == 1 ){
			whereto << " # Higgs potential" << endl;
		}
		
		whereto << "PoissSolnMethod = " << params->PoissSolnMethod << endl;
		if( params ->PoissSolnMethod == 0 ){
			whereto << " # not solving Poisson's equation " << endl;
		} 
		if( params ->PoissSolnMethod != 0 ){	

			if( params ->PoissSolnMethod == 1 ){
				whereto << " # FFT " << endl;
			} 
			if( params ->PoissSolnMethod == 2 ){
				whereto << " # relaxation " << endl;
				whereto << "PoissSolnRelaxMethod = " << params->PoissSolnRelaxMethod << endl;
				if ( params->PoissSolnRelaxMethod == 1 ){
					whereto << " # Gauss-Seidel" << endl;
				}
				if ( params->PoissSolnRelaxMethod == 2 ){
					whereto << " # Successive over relaxation (SOR)" << endl;
				}
				whereto << "Relaxation accuracy = " << params->PoissAccuracy << endl; 
			}
			
			whereto << "How often to solve Poisson's equation = " << params->PossSolveFreq << endl;
			
		}
		
		whereto << endl;
		
    } // END if( ID == 0 ){}
    
    // Write this at the end of a run
    if( ID == 100 ){
    
    	double TimeElapsed = params->TotalRunTime;
    	string TimeUnit = "milliseconds";
    	
    	// The default time-unit is milliseconds.
    	// If its longer than a second, write in seconds,
    	// and if longer than a minute, write in minutes, etc.
    	if(TimeElapsed > 1E3){
    	
	    	TimeElapsed = TimeElapsed / 1.0E3;
	    	TimeUnit = "seconds";
	    	
	    	if(TimeElapsed > 60){
	    		TimeElapsed = TimeElapsed / 60.0;
	    		TimeUnit = "minutes";
	    		
	    		if(TimeElapsed > 60){
	    			TimeElapsed = TimeElapsed / 60.0;
	    			TimeUnit = "hours";
    			} // END hour check
    			
    		} // END minute check
    		
    	} // END second check
    	
    	
	    whereto << endl;
    	if(params->flag == 0) whereto << "Completion without flag" << endl;
    	whereto << "Total run-time = " << TimeElapsed << " " << TimeUnit << endl;
    	whereto << endl;
    	whereto << "END" << endl;
    	whereto << endl;
    	
    } // END if( ID == 1 ){}
    
} // END PrintParams()




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF