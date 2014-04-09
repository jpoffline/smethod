
// timehistorystruct.h

#ifndef STRUCTTHIST_H
#define STRUCTTHIST_H

struct THIST{

	int timestep;
	double time;
	double val1;
	double val2;
	ofstream writeout;
	
	void write( struct THIST *timehistory ){
		
		timehistory->writeout << timehistory->timestep << " " ;
		timehistory->writeout << timehistory->time << " " ;
		timehistory->writeout << timehistory->val1 << " " ;		
		timehistory->writeout << timehistory->val2 << " " ;		
		timehistory->writeout << endl;
		
	} // END write()
	
	
};

#endif