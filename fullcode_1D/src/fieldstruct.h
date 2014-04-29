
// fieldstruct.h

// This contains all the field-type stuff.
// - compute finite difference derivatives
// - compute potential & derivative of the potential
// - compute equation of motion
// - update value of the field
// - items to dump to file

#ifndef STRUCTFIELD_H
#define STRUCTFIELD_H

#include "lapstencil.h"
#include "cosmologystruct.h"



struct FIELDCONTAINER{

	int ncom;
	double *vals, *deriv_x, *laplacian, *eom, *dpot, pot;
	double *FFTlap;
	
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
		
		
	// Get V(phi)
	void Getpot(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		// Get field at this location for easier reading of the code
		double *fld = new double[field->ncom];	
		for(int com = 0; com < field->ncom; com++){
		
			fld[com]=field->vals[ field->ind(grid->now,com,grid->loc_i,grid,field) ];
			
		}
				
		// Now construct the potentials
		// pottype can be used to select different potentials
		// ALSO must correspond to same pottype expressions in Getdpot() routine
		
		if(params->pottype == 0){
			// Massive scalar
			// V = m^2 phi^2 /2
			// ( note: potparam1 = m^2 )
			field->pot = 0.0;
			for(int com=0; com < field->ncom; com++){
			
				field->pot+= 0.5 * params->potparam1 * pow( fld[com] , 2.0 );
				
			}
			
		} // END pottype == 0
		
	
		if(params->pottype == 1){
			// Higgs potential
			// V = (phi^2 - 1)^2 / 4
			
			double ModPhiSq = 0.0;	// this will store |phi|^2
			field->pot = 0.0;
			// First, compute |phi|^2
			for(int com = 0; com < field->ncom; com++){
			
				ModPhiSq+= pow( fld[com] ,2.0);
				
			}
			
			field->pot = 0.25 * pow( ModPhiSq - 1.0 , 2.0 );
							
			
		} // END pottype == 1
		delete fld;
		
	} // END Getdpot()	
		
		
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////			
		
		
	// Get dV/dphi	
	void Getdpot(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){

		// Get field at this location for easier reading of the code		
		double *fld = new double[field->ncom];
		for(int com = 0; com < field->ncom; com++){
		
			fld[com]=field->vals[ field->ind(grid->now,com,grid->loc_i,grid,field) ];
			
		}
				
		// Now construct the derivative of the potential w.r.t each component of the field
		// pottype can be used to select different potentials
		// ALSO must correspond to same pottype expressions in Getpot() routine		
				
		
		if( params->pottype == 0 ){
			// Massive scalar
			// V = m^2 phi^2 /2
			
			for(int com = 0; com < field->ncom; com++){
			
				field->dpot[com] = params->potparam1 * fld[com];
				
			}
			
		} // END pottype == 0
		
		
		if( params->pottype == 1 ){
			// Higgs potential
			// V = (phi^2 - 1)^2 / 4
			double mod = 0.0; // this will store |phi|^2
			
			for(int com = 0; com < field->ncom; com++){
				mod+= pow( fld[com] ,2.0);
			}
			for(int com = 0; com < field->ncom; com++){
				field->dpot[com] = fld[com] * ( mod - 1.0 );
			}
			
		} // END pottype == 1
		
		/*
		// Derivative of scalar potential for Schrodinger-Poisson system.
		if( params->pottype == 2 ){
			
			// get scale factor
			double a = cosmology->a;
			
			// get solution to Poissons equation
			double V = poiss->V[grid->loc_i];
			
			// get hbar
			double hbar = cosmology->hbar;
			
			dpot[0] = pow( a , 2.0 ) * V / hbar * fld[1];
			dpot[1] = - pow( a , 2.0 ) * V / hbar * fld[0];		
		
		} // END params->pottype == 2
		*/
		
		// Deallocate memory
		delete fld;
		
	} // END Getdpot()
		
		
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////		
		
		
	// Get equation of motion	
	void GetEoM(struct DATA *params, struct FIELDCONTAINER *field, struct COSM *cosmology){

	// Returns the array holding the "E" parts
	// at this gridpoint -- returns for all components of the field.
	
		if( params->eomtype == 0  ){
			// Wave equation type:
			// E_i = nabla^2phi_i - dV/dphi_i, where i is the field component index
			for(int com = 0; com < field->ncom; com++){

				field->eom[com] = field->laplacian[com] - field->dpot[com];

			}
			
		}
		
		if( params->eomtype == 1 ){
			// Schrodinger type:
			// E_1 = - \nabla^2\phi_2 + dV/dphi_1
			// E_2 = \nabla^2\phi_1 - dV/dphi_2
			
			field->eom[0] = - cosmology->hbar / 2.0 * field->laplacian[1] + field->dpot[0];
			field->eom[1] = cosmology->hbar / 2.0 * field->laplacian[0] + field->dpot[1];			
			
	
			
		}
		
	} // END GetEoM()
	
	
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