
// relaxation.cpp

// This file contains the relaxation algorithms for solving Poisson's equation

#include "relaxation.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Gauss-Seidel algorithm

void gauss_seidel( struct POISS *poiss ){
    
    // Gauss-Seidel algorithm for solving 1D Poisson equation
    
    double h2 = poiss->h2;
    
    for(int i = 0; i < poiss->imax; i++){
    
        // Go & set ip & im, with periodic boundaries taken into account
    	poiss->SetrIP(i, poiss->imax, poiss);
		
		// Compute new value of V, via Gauss-Seidel
		poiss->rV[rVi(poiss->next, i)] = 0.5 * ( poiss->rV[ rVi(poiss->now, poiss->rip, poiss) ] 
										  + poiss->rV[ rVi(poiss->next, poiss->rim, poiss) ] 
										  - h2 * poiss->S[i] 
										);
            
    } // end i-loop
    
} // end gauss_seidel


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Successive over relaxation algorithm

void successive_over_relaxation( struct POISS *poiss){
    
    double h2 = poiss->h2;
    
    // SoR parameter
    double s = 2.0 / (1.0 + PI / poiss->imax );
    double halfs = 0.5 * s;
    
    for(int i = 0; i < poiss->imax; i++){
    
    	// Go & set ip & im, with periodic boundaries taken into account
	    poiss->SetrIP(i, poiss->imax, poiss);
	    
		// Compute new value of the potential V, via SOR
    	poiss->rV[rVi(poiss->next, i)] = ( 1.0 - s ) * poiss->rV[ rVi(poiss->now, i, poiss) ] 
    										+ halfs * ( poiss->rV[ rVi(poiss->now, poiss->rip, poiss) ] 
    													  + poiss->rV[ rVi(poiss->next, poiss->rim, poiss) ] 
    												   	  - h2 * poiss->S[i] 
    													);    
    			
    } // end i-loop
    
} // end successive_over_relaxation


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to return error on the relaxed solution.
// When this is smaller than desired accuracy, 
// the relaxation method stops iterating.

double relaxerror(struct POISS *poiss){
    
    // Relative value of the error
    double error = 0.0;
    
    // Number of points with non-zero value of the error
    int n = 0;
    
    double newpot, oldpot;
    
    for(int i = 0; i < poiss->imax; i++){
       
		newpot = poiss->rV[ rVi(poiss->next, i, params) ];
		
		if( newpot != 0 ){
			
			oldpot = poiss->rV[ rVi(poiss->now, i, params) ];
			
			if( newpot != oldpot ){
			
				error+= abs( 1.0 - newpot / oldpot );
				++n;
			
			} // end if()
			
		}// end if()
       
    }// end i-loop
    
    if( n != 0 ) error = error/n;
    
    return error;
    
} // end SolvePossion_error()


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF