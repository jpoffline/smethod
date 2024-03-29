
// gridstruct.h

// This holds everthing to do with the grid
// (including a copy of h & ht)
// Mainly used to hold the current time-step identifier
//  and + & - lattice site numbers, taking into account periodic boundaries

#ifndef STRUCTGRID_H
#define STRUCTGRID_H


struct GRIDINFO{

	int imin,imax;
	int loc_i;
	int ip,im;
	int t, prev,now,next;
	double h,h2,hh,htht,ht;
	int ntimsteps;

	// Routine to get time index correct
	// The algorithm leapfrogs through time
	void SetTime(int t, struct GRIDINFO *grid){
	
		grid->prev = 0;
		grid->now = 1;
		grid->t = t;
		if( t % 2 == 0 ){
		
			grid->prev = 1;
			grid->now = 0;
			
		}
		
		grid->next = grid->prev;
		
	} // END SetTime()


	// Routine to hold current lattice site
	//  and forwards/backwards lattice sites
	void GetPos(int pos, struct GRIDINFO *grid, int coord){
	
		// In the 3D code, "coord" corresponds to each of the 3 spatial
		// coordinates; here, only coord = 0 is used (that corresponds to "x")
		
		if(coord==0){
			
			// Store the current position, i
			grid->loc_i = pos;
			// Store ip = i + 1 (taking into account periodic boundaries)
			grid->ip = grid->GetP(pos,grid->imax);
			// Store im = i - 1 (taking in account periodic boundaries)
			grid->im = grid->GetM(pos,grid->imax);	
			
		}
		
	} // END GetPos()
	
	// Routine to return PLUS site: XP = X + 1
	//  Takes into account periodic boundary conditions
	int GetP(int pos, int posmax){
	
		int ret = pos + 1;
		if( ret < 0 ) ret = posmax - 1;
		if( ret > posmax - 1 ) ret = 0;
		return ret;
		
	} // END GetP()
	
	// Routine to return MINUS site: XM = X - 1
	//  Takes into account periodic boundary conditions
	int GetM(int pos, int posmax){
	
		int ret = pos - 1;
		if( ret < 0 ) ret = posmax - 1;
		if( ret > posmax - 1 ) ret = 0;
		return ret;
		
	} // END GetM()
	

};


#endif



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF