#ifndef STRUCTGRID_H
#define STRUCTGRID_H


struct GRIDINFO{

	int imin,imax,jmin,jmax,kmin,kmax;
	int loc_i,loc_j,loc_k;
	int ip,im,jp,jm,kp,km;
	int prev,now,next;
	double h,h2,hh,htht,ht;


	// Routine to get time index correct
	// The algorithm leapfrogs through time
	void SetTime(int t, struct GRIDINFO *grid){
	
		grid->prev = 0;
		grid->now = 1;
		
		if( t%2 == 0 ){
		
			grid->prev = 1;
			grid->now = 0;
			
		}
		
		grid->next = grid->prev;
		
	} // END SetTime()


	// Routine to hold current lattice site
	//  and forwards/backwards lattice sites
	void GetPos(int pos, struct GRIDINFO *grid, int coord){
	
		// Treats different coordinates differently (only due to different Xmax's
		//  (there is probably a better way of doing this)
		
		if(coord==0){
			grid->loc_i = pos;
			grid->ip = grid->GetP(pos,grid->imax);
			grid->im = grid->GetM(pos,grid->imax);	
		}
		if(coord==1){
			grid->loc_j = pos;
			grid->jp = grid->GetP(pos,grid->jmax);
			grid->jm = grid->GetM(pos,grid->jmax);	
		}
		if(coord==2){
			grid->loc_k = pos;
			grid->kp = grid->GetP(pos,grid->kmax);
			grid->km = grid->GetM(pos,grid->kmax);	
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