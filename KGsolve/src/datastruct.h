#ifndef STRUCTDATA_H
#define STRUCTDATA_H

struct DATA{

    double  h, ht;
    int     ntimsteps, imax, jmax, kmax;
    int     cmax;
    double  accuracy;
    int     derivsaccuracy; 
    int 	pottype, inittype, evoltype;
    int 	screenfreq, filefreq, thistfreq;
    double 	potparam1;
    int 	flag;
    
    int fx_i, fx_j, fx_k;
    
    double 	TotalRunTime;
    string 	OutDir, RunID;
    
}; // END DATA{}

#endif