
// timehistorystruct.h

#ifndef STRUCTTHIST_H
#define STRUCTTHIST_H
#include "poissstruct.h"

struct THIST{
	
	// Items to be stored as "time-history"
	int *timestep;
	double *time;
	double *poisserr;
	double *scalefactor;
	double *Hubble;
	double *eta;
	
	int thistfreq;

	
	ofstream TimeHistoryFile;

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	// Function to setup timehistory arrays.
	// The timehistory file is only written every
	// "thistfreq" time-steps, and everything inbetween is
	// saved in an array: this creates the right sized arrays.
	
	void setup( struct THIST *timehistory, int thistfreq){
	
		timehistory->thistfreq = thistfreq;
		timehistory->timestep = new int[thistfreq+1];
		timehistory->time  = new double[thistfreq+1];		
		timehistory->poisserr = new double[thistfreq+1];
		timehistory->scalefactor = new double[thistfreq+1];
		timehistory->Hubble = new double[thistfreq+1];	
		timehistory->eta = new double[thistfreq+1];			
			
	} // END setup()
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to set the items of the time-history.
	
	void SetItems(int th, struct GRIDINFO *grid, struct FIELDCONTAINER *field, struct THIST *timehistory){
			
		// Set the items to go into the timehistory
		
		timehistory->timestep[th] = grid->t;
		timehistory->time[th] = grid->t * grid->ht;
		if( field->poiss.method > 0 ) timehistory->poisserr[th] = field->poiss.poisserr;
		timehistory->eta[th] = field->cosmology.eta;
		timehistory->scalefactor[th] = field->cosmology.a;
		timehistory->Hubble[th] = field->cosmology.H;

			
	} // END SetTimeHistoryItems()
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	// Function to dump the time-history items to desired location.
	// This is the order that the items get dumped into file/screen.
	
	void TimeHistoryDump( ostream& whereto, struct THIST *timehistory, int th ){
		
		whereto << timehistory->timestep[th] << " " ;
		whereto << timehistory->time[th] << " " ;
		whereto << timehistory->poisserr[th] << " " ;
		whereto << timehistory->eta[th] << " " ;
		whereto << timehistory->scalefactor[th] << " " ;
		whereto << timehistory->Hubble[th] << " " ;
		
		// Make sure to write a new line after writing all the info you want.
		whereto << endl;
		
	} // END TimeHistoryDump()
	

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	// Function to write timehistory info to desired location.
	// This loops over all the "stored" values.
	
	void write( ostream& whereto, struct THIST *timehistory, int HowMany ){
		
		for(int th = 0; th <= HowMany; th++) TimeHistoryDump(whereto, timehistory, th);
		
	} // END write()
	
	
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	
	// Function to deallocate memory 
	
	void CleanUp( struct THIST *timehistory ){
	
		delete timehistory->timestep;
		delete timehistory->time;
		delete timehistory->poisserr;
		delete timehistory->scalefactor;
		delete timehistory->Hubble;		
		delete timehistory->eta;		
		
	} // END CleanUp()
	
};

#endif



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF