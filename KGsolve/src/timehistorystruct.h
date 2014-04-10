
// timehistorystruct.h

#ifndef STRUCTTHIST_H
#define STRUCTTHIST_H

struct THIST{

	int timestep,nFieldVals_thist;
	double time;

	ofstream writeout;
	
	// Arrays holding oordinates of the field
	// to be dumped into timehistory file
	int *val_i, *val_j, *val_k;
	double *val;	

	
	void SetFieldValDump( struct THIST *timehistory ) {
	
		// How many field values do you want to dump into time-history?
		timehistory->nFieldVals_thist = 2;
	
		timehistory->val = new double[timehistory->nFieldVals_thist];
		timehistory->val_i = new int[timehistory->nFieldVals_thist];
		timehistory->val_j = new int[timehistory->nFieldVals_thist];
		timehistory->val_k = new int[timehistory->nFieldVals_thist];	
		int vv = 0;	
		
		
		// Where are their coordinates?
		
		timehistory->val_i[vv] = 1;
		timehistory->val_j[vv] = 1;
		timehistory->val_k[vv] = 1;	
		vv++;
		timehistory->val_i[vv] = 1;
		timehistory->val_j[vv] = 5;
		timehistory->val_k[vv] = 1;	
	
	
	} // END SetFieldValDump()
	
		
	void write( struct THIST *timehistory ){
		
		// Write time-step number
		timehistory->writeout << timehistory->timestep << " " ;
		// Write physical time
		timehistory->writeout << timehistory->time << " " ;
		
		// Dump field values at the specified coordinates
		for(int v = 0; v< timehistory->nFieldVals_thist; v++){

			timehistory->writeout << timehistory->val[v] << " " ;	

		}	

		timehistory->writeout << endl;
		
	} // END write()
	
	void CleanUp( struct THIST *timehistory ){
	
		delete timehistory->val_i;
		delete timehistory->val_j;
		delete timehistory->val_k;
		delete timehistory->val;
	
	}
	
};

#endif