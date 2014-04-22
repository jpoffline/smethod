
// datastruct.h

// This contains the struct which holds all "parameters" (called as params in the main code)
// e.g. space & time step-sizes, number of time-steps

#ifndef STRUCTDATA_H
#define STRUCTDATA_H

struct DATA{

    double  h, ht;
    int     ntimsteps, imax;
    int     cmax;
    double  accuracy;
    int     derivsaccuracy; 
    int 	pottype, inittype, evoltype, eomtype;
    int 	field_lap_type;
    int 	screenfreq, filefreq, thistfreq;
    double 	potparam1;
    double 	initparam1;
    int 	flag;
    
    int lapstencilchoice;
    
    int PoissSolnMethod;
	int PoissSolnRelaxMethod;
	int PoissSourceType;
	double PoissAccuracy;
    
    double 	TotalRunTime;
    string 	OutDir, RunID;
    
}; // END DATA{}

#endif



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF