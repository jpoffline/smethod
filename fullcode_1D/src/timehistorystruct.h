
// timehistorystruct.h

#ifndef STRUCTTHIST_H
#define STRUCTTHIST_H

struct THIST{

	int *timestep;
	double *time;
	double *poisserr;
	int thistfreq;
	ofstream writeout;

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	// Function to setup timehistory arrays.
	// The timehistory file is only written every
	// "thistfreq" time-steps, and everything inbetween is
	// saved in an array
	
	void setup( struct THIST *timehistory, int thistfreq){
	
		timehistory->thistfreq = thistfreq;
		timehistory->timestep = new int[thistfreq+1];
		timehistory->time  = new double[thistfreq+1];		
		timehistory->poisserr = new double[thistfreq+1];
		
	} // END setup()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	// Function to write timehistory info to file
		
	void write( struct THIST *timehistory, int HowMany ){
		
		for(int th = 0; th <= HowMany; th++){
			timehistory->writeout << timehistory->timestep[th] << " " ;
			timehistory->writeout << timehistory->time[th] << " " ;
			timehistory->writeout << timehistory->poisserr[th] << " " ;
			timehistory->writeout << endl;
		}
		
	} // END write()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to deallocate memory 
	
	void CleanUp( struct THIST *timehistory ){
	
		delete timehistory->timestep;
		delete timehistory->time;
		delete timehistory->poisserr;
		
	} // END CleanUp()
	
};

#endif



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF