
// relaxation.cpp

// This file contains the relaxation algorithms for solving Poisson's equation

#include "relaxation.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Gauss-Seidel algorithm

void gauss_seidel( struct FIELDCONTAINER *field ){
    
    // Gauss-Seidel algorithm for solving 1D Poisson equation
    
    double h2 = field->poiss.h2;
    
    for(int i = 0; i < field->poiss.imax; i++){
    
        // Go & set ip & im, with periodic boundaries taken into account
    	field->poiss.SetrIP(i, field->poiss.imax, field);
		
		// Compute new value of V, via Gauss-Seidel
		field->poiss.rV[field->poiss.rVi(field->poiss.next, i, field)] 
										= 0.5 * ( field->poiss.rV[ field->poiss.rVi(field->poiss.now, field->poiss.rip, field) ] 
										  + field->poiss.rV[ field->poiss.rVi(field->poiss.now, field->poiss.rim, field) ] 
										  - h2 * field->poiss.S[i] 
										);
        // Using the "next" value for im needs to be done more carefully near the boundary!    
            
    } // end i-loop
    
} // end gauss_seidel


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Successive over relaxation algorithm

void successive_over_relaxation( struct FIELDCONTAINER *field ){
    
    double h2 = field->poiss.h2;
    
    // SoR parameter
    double s = 2.0 / (1.0 + PI / field->poiss.imax );
    double halfs = 0.5 * s;
    
    for(int i = 0; i < field->poiss.imax; i++){
    
    	// Go & set ip & im, with periodic boundaries taken into account
	    field->poiss.SetrIP(i, field->poiss.imax, field);
	    
		// Compute new value of the potential V, via SOR
    	field->poiss.rV[field->poiss.rVi(field->poiss.next, i, field)] 
										= ( 1.0 - s ) * field->poiss.rV[ field->poiss.rVi(field->poiss.now, i, field) ] 
    										+ halfs * ( field->poiss.rV[ field->poiss.rVi(field->poiss.now, field->poiss.rip, field) ] 
    													  + field->poiss.rV[ field->poiss.rVi(field->poiss.next, field->poiss.rim, field) ] 
    												   	  - h2 * field->poiss.S[i] 
    													);    
    			
    } // end i-loop
    
} // end successive_over_relaxation


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Function to return error on the relaxed solution.
// When this is smaller than desired accuracy, 
// the relaxation method stops iterating.

double relaxerror( struct FIELDCONTAINER *field ){
    
    // Relative value of the error
    double error = 0.0;
    
    // Number of points with non-zero value of the error
    int n = 0;
    
    double newpot, oldpot;
    
    for(int i = 0; i < field->poiss.imax; i++){
       
		newpot = field->poiss.rV[ field->poiss.rVi(field->poiss.next, i, field) ];
		
		if( newpot != 0 ){

			oldpot = field->poiss.rV[ field->poiss.rVi(field->poiss.now, i, field) ];
			
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