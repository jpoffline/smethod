
// fieldstruct.h

// This contains all the field-type stuff.
// - compute finite difference derivatives
// - compute potential & derivative of the potential
// - compute equation of motion
// - update value of the field
// - items to dump to file

#ifndef STRUCTFIELD_H
#define STRUCTFIELD_H

double DofPot(double field, struct DATA *params);

struct FIELDCONTAINER{

	int ncom;
	double *vals, *deriv_x, *deriv_y, *deriv_z, *laplacian, *eom, *dpot, pot;
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to return array index corresponding to 
	// time, component, and	spatial location.
	int ind(int t,int com,int i,int j,int k,struct GRIDINFO *grid,struct FIELDCONTAINER *field){
	
		return t * field->ncom * grid->imax * grid->jmax * grid->kmax
			   + com * grid->imax * grid->jmax * grid->kmax
			   + i * grid->jmax * grid->kmax
			   + j * grid->kmax
			   + k;	
		
	} // END ind()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to compute spatial derivatives
	// Choose 2nd or 4th order accurate
	void GetDeriv(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		if( params->derivsaccuracy == 2 ){
		
			field->GetDeriv_2(grid,field);
			
		}
		
		if( params->derivsaccuracy == 4 ){
		
			field->GetDeriv_4(grid,field);
			
		}
		
	} // END GetDeriv()
	
				
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////		
	
				
	// Second order accurate finite-difference spatial derivatives	
	void GetDeriv_2(struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		// df/dx = ( f[x+1] - f[x-1] ) / ( 2*h )
		// d2f/dx2 = ( f[x+1] + f[x-1] - 2 f[x] ) / ( h*h )
		// nabla^2f = d2f/dx2 + d2f/dy2 + d2f/dz2
	
		double f0,fip,fim,fjp,fjm,fkp,fkm;
		
		// Get derivatives for each component
		for(int com = 0;com < field->ncom; com++){
		
			// Get current value of the field
			f0 = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
			// Get phi(i+1)
			fip = field->vals[ field->ind(grid->now,com,grid->ip,grid->loc_j,grid->loc_k,grid,field) ];
			// Get phi(i-1)
			fim = field->vals[ field->ind(grid->now,com,grid->im,grid->loc_j,grid->loc_k,grid,field) ];
			// Get phi(j+1)
			fjp = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->jp,grid->loc_k,grid,field) ];		
			// Get phi(j-1)
			fjm = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->jm,grid->loc_k,grid,field) ];		
			// Get phi(k+1)
			fkp = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->kp,grid,field) ];				
			// Get phi(k-1)
			fkm = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->km,grid,field) ];	
	
			// dphi/fx			
			field->deriv_x[com] = ( fip - fim ) / grid->h2;
			// dphi/fy
			field->deriv_y[com] = ( fjp - fjm ) / grid->h2;
			// dphi/dz
			field->deriv_z[com] = ( fkp - fkm ) / grid->h2;
			// nabla^2 phi (3D)
			field->laplacian[com] = ( fip + fim + fjp + fjm + fkp + fkm - 6.0 * f0 ) / grid->hh;

		} // END com-loop

	} // END GetDeriv()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Fourth order accurate spatial derivatives
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
		
			fld[com]=field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
			
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
		
			fld[com]=field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
			
		}
				
		// Now construct the derivative of the potential w.r.t each component of the field
		// pottype can be used to select different potentials
		// ALSO must correspond to same pottype expressions in Getpot() routine		
				
		
		if(params->pottype == 0){
			// Massive scalar
			// V = m^2 phi^2 /2
			
			for(int com = 0; com < field->ncom; com++){
			
				field->dpot[com] = params->potparam1 * fld[com];
				
			}
			
		} // END pottype == 0
		
		
		if(params->pottype == 1){
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
		
		delete fld;
		
	} // END Getdpot()
		
		
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////		
		
		
	// Get equation of motion	
	void GetEoM(struct DATA *params, struct FIELDCONTAINER *field){

	// Returns the array holding the "E" parts
	// at this gridpoint -- returns for all components of the field.
	
		if( params->eomtype == 0  ){
			// Wave equation type:
			// E = nabla^2phi - dV/dphi
			for(int com = 0; com < field->ncom; com++){

				field->eom[com] = field->laplacian[com] - field->dpot[com];

			}
			
		}
		if( params->eomtype == 1 ){
			// Schrodinger type:

			
		}
		
	} // END GetEoM()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to update field value from 2nd order EoM
	void UpdateField(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		double fp,fn;
		
		for(int com = 0; com < field->ncom; com++){
		
			// Get previous value of the field
			fp = field->vals[ field->ind(grid->prev,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
		
			// Get current value of the field
			fn = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
		
			// Update value of the field: choose which rule to use via evoltype
			
			if( params->evoltype == 0 ){
			
				// (1) Gradient flow
			
				fp = field->eom[com] * grid->ht + fn;
			
			}
			
			if( params->evoltype == 1 ){
			
				// (2) 2nd order Klein-Gordon
			
				fp = field->eom[com] * grid->htht - fp + 2.0 * fn;
			
			}
			
			// Dump computed field into the "new" value of the field
			
			field->vals[ field->ind(grid->next,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ] = fp;
			
		}
		
	} // END UpdateField()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	void WriteFieldData(ostream& whereto, struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		// At the moment, this is still outputting quite a lot more information than a proper run would need
		whereto <<  grid->loc_i << " " <<  grid->loc_j << " " <<  grid->loc_k << " " ;
		whereto <<  grid->ip << " " <<  grid->jp << " " <<  grid->kp << " " ;
		whereto <<  grid->im << " " <<  grid->jm << " " <<  grid->km << " " ;
	
		for(int com=0;com < field->ncom; com++){
		
			// Output current value of the field
			whereto << field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ] << " " ;
			
			// Output previous value of the field
			whereto << field->vals[ field->ind(grid->prev,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ] << " " ;			
			
		}
		// Print newline
		whereto << endl;
		
	} // END WriteFieldData()
	
	
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////	
	
	
	// Routine to delete any arrays that were allocated
	void CleanField(struct FIELDCONTAINER *field){
	
		delete field->vals;		
		delete field->laplacian;
		delete field->deriv_x;
		delete field->deriv_y;
		delete field->deriv_z;
		delete field->eom;
		delete field->dpot;	 
		
	} // END CleanField()
	
}; // END FIELDCONTAINER{}

#endif