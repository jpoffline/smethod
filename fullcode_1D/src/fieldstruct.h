
// fieldstruct.h

// This contains all the routines for field-type stuff.
// - compute finite difference derivatives
// - update value of the field
// - items to dump to file

#ifndef STRUCTFIELD_H
#define STRUCTFIELD_H

struct FIELDCONTAINER{
	
	
	
	int ncom;
	double *vals, *deriv_x, *laplacian, *eom, *dpot, pot;
	double *FFTlap;
	
	// Two-subclasses:
	// (1) containing all background cosmology
	#include "cosmologystruct.h"	
	COSM cosmology;
	
	// (2) containing all fields and routines for solving
	// Poissons equation
	#include "poissstruct.h"
	POISS poiss;
	
	// (3) containing all the interesting physics equations
	#include "eqnsstruct.h"
	EQUATIONS equations;
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to return array index corresponding to 
	// time, component, and	spatial location.
	int ind(int t,int com,int i,struct GRIDINFO *grid,struct FIELDCONTAINER *field){
	
		return t * field->ncom * grid->imax
			   + com * grid->imax
			   + i; 
		
	} // END ind()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to compute spatial derivatives
	// Choose FFT or second order finite difference,
	// and then, if finite difference, 2nd or 4th order accurate
	void GetDeriv(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		// Finite difference derivatives
		if( params->field_lap_type == 0 ){
			
			// Second order accurate
			if( params->derivsaccuracy == 2 ) 
				field->GetDeriv_2(grid,field, params);
			
			// Fourth order accurate
			if( params->derivsaccuracy == 4 ) 
				field->GetDeriv_4(grid,field);
			
		}
		
		// FFT derivatives
		if( params->field_lap_type == 1 )
			for(int c = 0; c < 	field->ncom; c++)
				field->laplacian[c] = field->FFTlap[field->ind(grid->now,c,grid->loc_i,grid,field)];
		
	} // END GetDeriv()
	
				
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////		
	
				
	// Second order accurate finite-difference spatial derivatives	
	void GetDeriv_2(struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct DATA *params){
	
		// df/dx = ( f[x+1] - f[x-1] ) / ( 2*h )
		// d2f/dx2 = ( f[x+1] + f[x-1] - 2 f[x] ) / ( h*h )
		// nabla^2f = d2f/dx2 
	
		double f0, fip, fim;
		
		// Get derivatives for each component
		for(int com = 0; com < field->ncom; com++){
		
			f0 = field->vals[ field->ind(grid->now,com,grid->loc_i,grid,field) ];
			fip = field->vals[ field->ind(grid->now,com,grid->ip,grid,field) ];
			fim = field->vals[ field->ind(grid->now,com,grid->im,grid,field) ];

			// dphi/fx			
			field->deriv_x[com] = ( fip - fim ) / grid->h2;

			// nabla^2 phi 
			if(params->field_lap_type == 0){
			
				field->laplacian[com] = ( fip + fim - 2.0 * f0 ) / grid->hh;
			
			}
			
			
			
		} // END com-loop

	} // END GetDeriv()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Fourth order accurate finite-difference spatial derivatives
	void GetDeriv_4(struct GRIDINFO *grid, struct FIELDCONTAINER *field){
		
		// This can easily be coded up
		// BUT: would require an additional rethink for the GetM etc routines
		//	since more gridpoints are required to compute the derivatives.
		
	} // END GetDeriv_4()
	
		
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
			
	
	// Routine to update field value from 2nd order EoM, or gradient flow
	void UpdateField(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		double fp,fn;
		
		for(int com = 0; com < field->ncom; com++){
		
			// Get previous value of the field
			fp = field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
		
			// Get current value of the field
			fn = field->vals[ field->ind(grid->now,com,grid->loc_i,grid,field) ];
		
			// Update value of the field: choose which rule to use via evoltype
			
			// (1) Gradient flow
			if( params->evoltype == 0 ){
			
				fp = field->eom[com] * grid->ht + fn;
			
			}
			
			// (2) 2nd order wave equation
			if( params->evoltype == 1 ){				
			
				fp = field->eom[com] * grid->htht - fp + 2.0 * fn;
			
			}
			
			// Dump computed field into the "new" value of the field
			
			field->vals[ field->ind(grid->next,com,grid->loc_i,grid,field) ] = fp;
			
		}
		
	} // END UpdateField()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	void WriteFieldData(ostream& whereto, struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		// At the moment, this is still outputting quite a lot more information than a proper run would need
		whereto <<  params->h * grid->loc_i <<  " " ;
	
		for(int com = 0;com < field->ncom; com++){
			// Output current value of the components of the field
			whereto << field->vals[ field->ind(grid->now,com,grid->loc_i,grid,field) ] << " " ;
		}
		// Print newline
		whereto << endl;
		
	} // END WriteFieldData()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to delete any arrays that were allocated
	void CleanField(struct FIELDCONTAINER *field){
	
		delete field->vals;	
		delete field->FFTlap;	
		delete field->laplacian;
		delete field->deriv_x;
		delete field->eom;
		delete field->dpot;	 
		
	} // END CleanField()
	
}; // END FIELDCONTAINER{}

#endif



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF