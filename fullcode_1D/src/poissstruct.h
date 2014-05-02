
// poissstruct.h


#ifndef STRUCTPOISS_H
#define STRUCTPOISS_H

struct POISS{

	// nabla^2 V = S

	int imax;	
	int now, next;
	int rip, rim;
	double h2, h;
	
	// ID of method to solve Poisson's equation
	// 1 = FFT
	// 2 = relaxation
	int method;
	
	// ID to specify which relaxation method to use
	// 1 = Gauss-Siedel
	// 2 = SOR
	int relaxmethod;
	double poisserr;
	
	// ID to specify the type of source of the Poisson equation
	int source_type;
	
	// Accuracy with which to solve Poisson equation
	double accuracy;
	
	// Array to hold source of Poisson equation
    double *S;
    // Array to hold potential to be solved for
    double *V;
    
    // Array to hold potential for relaxation algorithms
    double *rV;
    
	// How often do we want to solve Poisson's equation?
	int PossSolveFreq;
    
    ////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	
    // Function to set the "now" and "next" for V
    // as used by the relaxation algorithms
    
    void SetStep( int step, struct FIELDCONTAINER *field ){
    
		field->poiss.now = 1;
		field->poiss.next = 0;
	
		if( step % 2 == 0 ){
		
			field->poiss.now = 0; 
			field->poiss.next = 1;
			
		}
		
    } // END SetStep()
	
	
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	
	// Function to return array index of V used for relaxation
	
	int rVi(int t, int i, struct FIELDCONTAINER *field ){
	
		return t * field->poiss.imax + i;
	
	} // END rVi()
	
	
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	
	// Function to set ip = i + 1 and im = i - 1 for relaxation
	
	void SetrIP(int i, int posmax, struct FIELDCONTAINER *field){

		int id;
		
		id = i + 1;
		if( id == posmax ) id = 0;
		field->poiss.rip = id;
		
		id = i - 1;
		if( id < 0 ) id = posmax - 1;
		field->poiss.rim = id;
		
	} // END SetrIP()
	
	
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
    
    // Function to setup the source of the Poisson equation
    
    void SetupSource(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
		
		// Note, \nabla^2V = S, and the form of S is set below
		
		if( field->poiss.source_type == 1 ){	
			
			// Set S = |psi|^2 - 1
			
			double modSq;
			for(int i = 0; i < field->poiss.imax; i++ ){
				modSq = 1.0;
				for(int com = 0; com < params->cmax; com++)
					modSq*=pow(field->vals[field->ind(grid->now,com,i,grid,field)],2.0);

				field->poiss.S[i] = modSq - 1.0;
			}
			
		} // END poiss->source_type == 1
		
		if( field->poiss.source_type == 2 ){
			
			// Set S = 3 H0 / ( 2a ) (|psi|^2 - 1)
			
			double modSq;
			for(int i = 0; i < field->poiss.imax; i++ ){
				modSq = 1.0;
				for(int com = 0; com < params->cmax; com++)
					modSq*=pow(field->vals[field->ind(grid->now,com,i,grid,field)],2.0);
				
				field->poiss.S[i] = 3.0 * field->cosmology.H0 / 2.0 / field->cosmology.a * ( modSq - 1.0 );
			}
			
		} // END poiss->source_type == 2

		
	} // END SetupSource()
    
    
    ////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
    
    // Function to write Poisson solution info to file
    
    void WritePoissData(int i, ostream& whereto, struct POISS *poiss ){
	
	    whereto << poiss->h * i << " " << poiss->V[ i ] << endl;
		
	} // END WriteFieldData()
    
    
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    
    // Function to deallocate arrays allocated for the Poisson equation
    
    void CleanPoiss(struct POISS *poiss){
    
    	delete S;
    	delete V;
    	delete rV;
    
    } // END CleanPoiss()
    
    
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    
}; // END POISS{}

#endif



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF