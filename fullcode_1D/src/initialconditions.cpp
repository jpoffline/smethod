
// initialconditions.cpp

#include "initialconditions.h"

void SetInitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct COSM *cosmology){

	// This routine contains the different types of initial conditions
	// This is called by InitialConditions() -- at the end of this file
	//  at each lattice site

	int pos;
	double sign, ran, value;

	// homogeneous
	if( params->inittype == 0 ){
		for(int com = 0; com < field->ncom; com++){
			
			pos = field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos] = 0.001;
			
			pos = field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
			
		}
	} // END if( params->inittype==1 )
	
	// random distribution about the origin
	if( params->inittype == 1 ){
	
	//	srand(time(NULL));
		for(int com = 0; com < field->ncom; com++){
		
			ran=rand()/(double)RAND_MAX;
			sign = 1.0;
			if(ran<0.5){sign=-1.0;}
			ran=sign*rand()/(double)RAND_MAX;
			
			pos=field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos]=ran;
			
			pos=field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
		}
		
	} // END if( params->inittype == 1 )
	
	
	// square kink
	if( params->inittype == 2 ){
		for(int com = 0; com < field->ncom; com++){
		
			value = 1.0;
			if(grid->loc_i < 0.5 * grid->imax) value = -1.0;
			
			pos=field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos]=value;
			
			pos=field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
			
		}
	} // END if( params->inittype == 2 )
	
	
	// sine wave
	if( params->inittype == 3){
	
		double omega = params->initparam1 * 2.0 * PI / (  params->h * params->imax );
		double x = grid->loc_i * params->h;	

	
		value = sin( x * omega );
		for(int com = 0; com < field->ncom; com++){
			pos=field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos]=value;
			
			pos=field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
		}
	
	} // END if( params->inittype == 3 )
	
	// collapse initial conditions
	if (params->inittype == 4){
		
		// ci = sqrt(-1):
		dcmplx ci(0.0,1.0);
		
		// parameters
		double x = grid->loc_i * params->h;
		
		// derived quantities
		// scale factor (already set at start of run)
		
		double a = cosmology->a;
		
		// overdensity
		double delta = a * cos( PI * x / cosmology->L ); 
		double n = 1.0 + delta;
		double Vd = - 3.0/2.0 * cosmology->H0 * pow( cosmology->L / PI , 2.0 ) / a * delta;
		double phi = - 2.0 / 3.0 / cosmology->getH( a , cosmology ) * Vd; 
		
		
		
		dcmplx psi = sqrt(n) * exp( ci * phi / cosmology->hbar );
		
		// get "real" part of scalar field
		double fld_real = real(psi);		
		field->vals[ field->ind(grid->prev,0,grid->loc_i,grid,field) ] = fld_real;
		// get "imaginary" part of scalar field
		double fld_imag = imag(psi);
		field->vals[ field->ind(grid->prev,1,grid->loc_i,grid,field) ] = fld_imag;


		
	} // END if (params->inittype == 4)
	
} // END SetInitialConditions()




void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct COSM *cosmology){
	
	grid->SetTime(0,grid);

	for(int i = grid->imin; i < grid->imax; i++){
		grid->GetPos(i,grid,0);	

		SetInitialConditions(params,grid,field,cosmology);

	} // END i-loop

} // END InitialConditions





////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF