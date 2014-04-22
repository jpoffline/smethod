
// initialconditions.cpp

#include "initialconditions.h"

void SetInitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){

	// This routine contains the different types of initial conditions
	// This is called by InitialConditions() -- at the end of this file
	//  at each lattice site

	int pos;
	double sign, ran, value;

	if( params->inittype == 0 ){
		for(int com=0; com < field->ncom; com++){
			
			pos = field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos] = 0.001;
			
			pos = field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
			
		}
	} // END if( params->inittype==1 )
	
	if( params->inittype == 1 ){
	
	//	srand(time(NULL));
		for(int com=0; com < field->ncom; com++){
		
			ran=rand()/(double)RAND_MAX;
			sign = 1.0;
			if(ran<0.5){sign=-1.0;}
			ran=sign*rand()/(double)RAND_MAX;
			
			pos=field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos]=ran;
			
			pos=field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
		}
		
	} // END if( params->inittype==1 )
	
	if( params->inittype == 2 ){
		for(int com=0; com < field->ncom; com++){
		
			value = 1.0;
			if(grid->loc_i < 0.5 * grid->imax) value = -1.0;
			
			pos=field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos]=value;
			
			pos=field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
			
		}
	} // END if( params->inittype==2 )
	
	
	if( params->inittype == 3){
	
		double omega = params->initparam1 * 2.0 * PI / (  params->h * params->imax );
		double x = grid->loc_i * params->h;	

	
		value = sin( x * omega );
		for(int com=0; com < field->ncom; com++){
			pos=field->ind(grid->prev,com,grid->loc_i,grid,field);
			field->vals[pos]=value;
			
			pos=field->ind(grid->now,com,grid->loc_i,grid,field);
			field->vals[pos]=field->vals[ field->ind(grid->prev,com,grid->loc_i,grid,field) ];
		}
	
	}
	
	
} // END SetInitialConditions()




void InitialConditions(struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field){
	
	grid->SetTime(0,grid);

	for(int i = grid->imin; i < grid->imax; i++){
		grid->GetPos(i,grid,0);	

		SetInitialConditions(params,grid,field);

	} // END i-loop

} // END InitialConditions





////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF