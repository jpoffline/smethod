
// fieldstruct.h

#ifndef STRUCTFIELD_H
#define STRUCTFIELD_H

double DofPot(double field, struct DATA *params);

struct FIELDCONTAINER{

	int ncom;
	double *vals;	
	double *deriv_x,*deriv_y,*deriv_z,*laplacian,*eom,*dpot,*pot;
	// Routine to return array index corresponding to 
	// time, component, and	spatial location.
	int ind(int t,int com,int i,int j,int k,struct GRIDINFO *grid,struct FIELDCONTAINER *field){
	
		int	dum;
		dum = t * field->ncom * grid->imax * grid->jmax * grid->kmax;
		dum+= com * grid->imax * grid->jmax * grid->kmax;
		dum+= i * grid->jmax * grid->kmax;
		dum+= j * grid->kmax;
		dum+= k;
		return dum;
		
	}
	
	// Routine to compute spatial derivatives
	// Choose 2nd or 4th order accurate
	void GetDeriv(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		if(params->derivsaccuracy==2){
			field->GetDeriv_2(grid,field);
		}
		if(params->derivsaccuracy==4){
			field->GetDeriv_4(grid,field);
		}
		
	} // END GetDeriv()
				
	// Second order accurate spatial derivatives	
	void GetDeriv_2(struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		double f0,fip,fim,fjp,fjm,fkp,fkm;
	
		for(int com=0;com < field->ncom; com++){
		
			f0 = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
			fip = field->vals[ field->ind(grid->now,com,grid->ip,grid->loc_j,grid->loc_k,grid,field) ];
			fim = field->vals[ field->ind(grid->now,com,grid->im,grid->loc_j,grid->loc_k,grid,field) ];
			fjp = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->jp,grid->loc_k,grid,field) ];		
			fjm = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->jm,grid->loc_k,grid,field) ];		
			fkp = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->kp,grid,field) ];				
			fkm = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->km,grid,field) ];				
			field->deriv_x[com] = ( fip - fim ) / grid->h2;
			field->deriv_y[com] = ( fjp - fjm ) / grid->h2;
			field->deriv_z[com] = ( fkp - fkm ) / grid->h2;
			field->laplacian[com] = ( fip + fim + fjp + fjm + fkp + fkm - 6.0 * f0 ) / grid->hh;

		}

	} // END GetDeriv()
	
	// Fourth order accurate spatial derivatives
	void GetDeriv_4(struct GRIDINFO *grid, struct FIELDCONTAINER *field){
		
	}	
		
		
	// Get dV/dphi	
	void Getpot(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){

		double mod = 0.0;
		
		double *fld = new double[field->ncom];
		
		// Get field at this location for easier reading of the code
		for(int com = 0; com < field->ncom; com++){
			fld[com]=field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
		}
				
		
		if(params->pottype == 0){
			// Massive scalar
			for(int com=0; com < field->ncom; com++){
				field->pot[com] = 0.5 * params->potparam1 * pow( fld[com] , 2.0 );
			}
		} // END pottype == 0
		
		if(params->pottype == 1){
			// Higgs potential
			for(int com = 0; com < field->ncom; com++){
				mod+= pow( fld[com] ,2.0);
			}
			for(int com=0;com< field->ncom; com++){
				field->pot[com] = 0.25 * pow( mod - 1.0 , 2.0 );
			}
		} // END pottype == 1
		delete fld;
		
	} // END Getdpot()	
		
	// Get dV/dphi	
	void Getdpot(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){

		double mod = 0.0;		
		double *fld = new double[field->ncom];
		
		// Get field at this location for easier reading of the code
		for(int com = 0; com < field->ncom; com++){
			fld[com]=field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
		}
				
		
		if(params->pottype == 0){
			// Massive scalar
			for(int com=0; com < field->ncom; com++){
				field->dpot[com] = params->potparam1 * fld[com];
			}
		} // END pottype == 0
		
		if(params->pottype == 1){
			// Higgs potential
			for(int com = 0; com < field->ncom; com++){
				mod+= pow( fld[com] ,2.0);
			}
			for(int com=0;com< field->ncom; com++){
				field->dpot[com] = fld[com] * ( mod - 1.0 );
			}
		} // END pottype == 1
		delete FldL;
		
	} // END Getdpot()
		
		
		
	// Get EoM	
	void GetEoM(struct FIELDCONTAINER *field){
	
		for(int com=0; com < field->ncom; com++){
			field->eom[com] = field->laplacian[com] - field->dpot[com];
		}
		
	} // END GetEoM()
	
	// Routine to update field value from 2nd order EoM
	void UpdateField(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
		double fp,fn;
		for(int com=0;com < field->ncom; com++){
			fp = field->vals[ field->ind(grid->prev,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
			fn = field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ];
			if( params->evoltype == 0 ){
				fp = field->eom[com] * grid->ht + fn;
			}
			if( params->evoltype == 1 ){
				fp = field->eom[com] * grid->htht - fp + 2.0 * fn;
			}
			field->vals[ field->ind(grid->next,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ] = fp;
		}
		
	} // END UpdateField()
	
	void WriteFieldData(ostream& whereto, struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
		whereto <<  grid->loc_i << " " <<  grid->loc_j << " " <<  grid->loc_k << " " ;
		whereto <<  grid->ip << " " <<  grid->jp << " " <<  grid->kp << " " ;
		whereto <<  grid->im << " " <<  grid->jm << " " <<  grid->km << " " ;
		for(int com=0;com < field->ncom; com++){
			whereto << field->vals[ field->ind(grid->now,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ] << " " ;
			whereto << field->vals[ field->ind(grid->prev,com,grid->loc_i,grid->loc_j,grid->loc_k,grid,field) ] << " " ;			
		}
		whereto << endl;

		
	} // END WriteFieldData()
	
	
	// Routine to delete any arrays that were allocated
	void CleanField(struct FIELDCONTAINER *field){
	
		delete field->vals;		
		delete field->laplacian;
		delete field->deriv_x;
		delete field->deriv_y;
		delete field->deriv_z;
		delete field->eom;
		delete field->dpot;
		delete field->pot;
		
	} // END CleanField()
	
};

#endif