
// eqnsstruct.h

// Structure to contain all the model specific equations
// e.g. specific equations of motion, potentials

#ifndef STRUCTEQNS_H
#define STRUCTEQNS_H

struct EQUATIONS{
	
	
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
			double a = field->cosmology.a;
			
			// get solution to Poissons equation
			double V = poiss->V[grid->loc_i];
			
			// get hbar
			double hbar = field->cosmology.hbar;
			
			dpot[0] = pow( a , 2.0 ) * V / hbar * fld[1];
			dpot[1] = - pow( a , 2.0 ) * V / hbar * fld[0];		
		
		} // END params->pottype == 2
		*/
		
		// Deallocate memory
		delete fld;
		
	} // END Getdpot()
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
		
	// Get equation of motion, E.
	// E = \dot{phi} or E = \ddot{\phi}, depending on evoltype setting.
	// This returns "E".
	void GetEoM(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){

	// Returns the array holding the "E" parts
	// at this gridpoint -- returns for all components of the field.
	
		if( params->eomtype == 0  ){
			// Wave equation type:
			// E_i = nabla^2phi_i - dV/dphi_i, where i is the field component index
			for(int com = 0; com < field->ncom; com++)
				field->eom[com] = field->laplacian[com] - field->dpot[com];
			
		}
		
		if( params->eomtype == 1 ){
			// Schrodinger type:
			// E_1 = - \nabla^2\phi_2 + dV/dphi_1
			// E_2 = \nabla^2\phi_1 - dV/dphi_2
			
			field->eom[0] = - field->cosmology.hbar / 2.0 * field->laplacian[1] + field->dpot[0];
			field->eom[1] = field->cosmology.hbar / 2.0 * field->laplacian[0] + field->dpot[1];			

		}
 		
		if( params->eomtype == 2 ){
			// wave equation type:
			// with Poisson solution, V.
			// E_i = (1+4V)\nabla^2phi_i - a^2(1 + 2V) dU/dphi_i
 			double NewtPot = field->poiss.S[grid->loc_i];

			for(int com = 0; com < field->ncom; com++)
				field->eom[com] = ( 1.0 + 4.0 * NewtPot ) * field->laplacian[com] 
								- pow( field->cosmology.a , 2.0 ) * ( 1.0 + 2.0 * NewtPot ) * field->dpot[com];			

		}
 
				
	} // END GetEoM()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
};

#endif