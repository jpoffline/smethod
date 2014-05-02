
// setup.cpp

// This contains routines to read in the params.ini file
//  and populate items in the structs

#include "setup.h"

void Setup( struct DATA *params, struct GRIDINFO *grid, struct FIELDCONTAINER *field ){

	// Main caller routine to set everything up

	SetupGrid(grid, params);
	SetupField(params, field);
	SetupPoisson(params, field);
	SetupCosmology(params, field);
	
} // END Setup()

void SetupGrid(struct GRIDINFO *grid, struct DATA *params){

	grid->imin = 0;
	grid->imax = params->imax;
		
	grid->h = params->h;
	grid->hh = pow(params->h,2.0);
	grid->h2 = 2.0 * params->h;
	grid->ht = params->ht;
	grid->htht = pow(params->ht,2.0);
	grid->ntimsteps = params->ntimsteps;
	
} // END SetupGrid()

void SetupField(struct DATA *params, struct FIELDCONTAINER *field){

	field->ncom = params->cmax;
	field->vals = new double[ 2 * field->ncom * params->imax ];
	if( params->field_lap_type == 1 ) field->FFTlap = new double[ 2 * field->ncom * params->imax ];
	field->deriv_x = new double[ params->cmax ];
	field->laplacian = new double[ params->cmax ];
	field->eom = new double[ params->cmax ];
	field->dpot = new double [ params->cmax ];
	
} // END SetupField()

void SetupPoisson(struct DATA *params, struct FIELDCONTAINER *field){

	field->poiss.imax = params->imax;
	field->poiss.h = params->h;
	field->poiss.h2 = pow(params->h,2.0);
	field->poiss.method = params->PoissSolnMethod;
	field->poiss.relaxmethod = params->PoissSolnRelaxMethod;
	field->poiss.source_type = params->PoissSourceType;
	field->poiss.accuracy = params->PoissAccuracy;
	field->poiss.PossSolveFreq = params->PossSolveFreq;
    field->poiss.S = new double[ field->poiss.imax ];
    field->poiss.V = new double[ field->poiss.imax ];
    if( field->poiss.method == 2) field->poiss.rV = new double[ 2 * field->poiss.imax ];

} // END SetupPoisson()

void SetupCosmology(struct DATA *params, struct FIELDCONTAINER *field){
	
	// How to compute "L" from the number of grid points,
	// and space step-size?
	// NOTE: physical size of the box is imax*h
	field->cosmology.L = 0.5 * params->imax * params->h;
	
	// Get the current value of the Hubble parameter from what 
	// was read in from the params.ini file
	field->cosmology.H0 = params->H0;
	
	// Get the value of "hbar" from what 
	// was read in from the params.ini file
	field->cosmology.hbar = params->hbar;

} // END SetupCosmology()

void GetParams(int argc, char* argv[], struct DATA *params){

  	// Get the parameter file:
    IniReader inifile;
    // If the user specified a param file at runtime
    //  then use that, else, use the default, params.ini
    //   (User specifies via ./EXE my_param.par)
    
    if (argc > 1)
        inifile.read(argv[1]);
    else
        inifile.read("params.ini");


	// The numbers shown are the "default" values,
	// if nothing can be found in the params.ini file.
	
    params->h = inifile.getiniDouble("h",0.1);
    params->ht = inifile.getiniDouble("ht",0.01);
    params->derivsaccuracy = int(inifile.getiniDouble("derivsaccuracy",2));
    
    params->imax = int(inifile.getiniDouble("imax",10));
    params->cmax = int(inifile.getiniDouble("cmax",2));
    params->ntimsteps = int(inifile.getiniDouble("ntimsteps",100));
  
	params->pottype = int(inifile.getiniDouble("pottype",1));
	params->inittype = int(inifile.getiniDouble("inittype",1));
	params->evoltype = int(inifile.getiniDouble("evoltype",0));
	params->field_lap_type = int(inifile.getiniDouble("field_lap_type",0));
	params->eomtype = int(inifile.getiniDouble("eomtype",0));
	
	params->screenfreq = int(inifile.getiniDouble("screenfreq",10));
	params->filefreq = int(inifile.getiniDouble("filefreq",10));
	params->thistfreq = int(inifile.getiniDouble("thistfreq",1));
	params->potparam1 = inifile.getiniDouble("potparam1",0.0);
	params->initparam1 = inifile.getiniDouble("initparam1",0.0);
	
	params->OutDir = inifile.getiniString("OutDir","output/");
	params->RunID = inifile.getiniString("RunID","run_01");
	
	params->PoissSolnMethod = int(inifile.getiniDouble("PoissSolnMethod",1));
	params->PoissSolnRelaxMethod = int(inifile.getiniDouble("PoissSolnRelaxMethod",1));
	params->PoissSourceType = int(inifile.getiniDouble("PoissSourceType",1));
	params->PoissAccuracy =  inifile.getiniDouble("PoissAccuracy",1);
	params->PossSolveFreq =  inifile.getiniDouble("PossSolveFreq",1);
	
	params->H0 = inifile.getiniDouble("H0",1.0);
	params->hbar = inifile.getiniDouble("hbar",1.0);
	
	
	// Make sure output directory exists.
	// If it dosnt, this function will create it & tell you
	// it got created
    checkdirexists(params->OutDir);
	
	params->flag = 0;
	
	
} // END GetParams()

void CheckParams(struct DATA *params){

	// This checks the parameters for sanity reasons
	// If anything has been setup pathologically, code will not run
	
	if( params->imax < 0 ) params->flag = 1;
	
	if( params->flag != 0 ){
		cout << endl;	
		cout << "BAD PARAMETERS!" << endl;
		cout << "flag = " << params->flag << endl;
		cout << "TERMINATING CODE" << endl;
		cout << endl;
		
	}
	
} // END CheckParams()

void checkdirexists(string dir){
    
    using namespace boost::filesystem;
    
    if (!exists(dir + "/")) {
        cout << endl;
        cout << " --> Creating output directory" << endl;
        cout << endl;
        
        create_directory(dir);
    }
    
} // END checkdirexists()



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// EOF