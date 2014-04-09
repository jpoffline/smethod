
// setup.cpp

#include "setup.h"

void SetupGrid(struct GRIDINFO *grid, struct DATA *params){

	grid->imin = 0;
	grid->jmin = 0;
	grid->kmin = 0;
	grid->imax = params->imax;
	grid->jmax = params->jmax;
	grid->kmax = params->kmax;
	
	grid->h = params->h;
	grid->hh = pow(params->h,2.0);
	grid->h2 = 2.0 * params->h;
	grid->ht = params->ht;
	grid->htht = pow(params->ht,2.0);
	
} // END SetupGrid()

void SetupField(struct DATA *params, struct FIELDCONTAINER *field){

	field->ncom = params->cmax;
	field->vals = new double[ 2 * field->ncom * params->imax * params->jmax * params->kmax ];
	field->deriv_x = new double[ params->cmax ];
	field->deriv_y = new double[ params->cmax ];
	field->deriv_z = new double[ params->cmax ];
	field->laplacian = new double[ params->cmax ];
	field->eom = new double[ params->cmax ];
	field->dpot = new double [ params->cmax ];
	field->pot = new double [ params->cmax ] ;
	
} // END SetupField()

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

    params->h = inifile.getiniDouble("h",0.1);
    params->ht = inifile.getiniDouble("ht",0.01);
    params->accuracy = inifile.getiniDouble("accuracy",0.0001);
    params->derivsaccuracy = int(inifile.getiniDouble("derivsaccuracy",2));
    
    params->imax = int(inifile.getiniDouble("imax",10));
    params->jmax = int(inifile.getiniDouble("jmax",10));
    params->kmax = int(inifile.getiniDouble("kmax",10));
    params->cmax = int(inifile.getiniDouble("cmax",2));
    params->ntimsteps = int(inifile.getiniDouble("ntimsteps",100));
    
	params->pottype = int(inifile.getiniDouble("pottype",1));
	params->inittype = int(inifile.getiniDouble("inittype",1));
	params->evoltype = int(inifile.getiniDouble("evoltype",0));
	params->screenfreq = int(inifile.getiniDouble("screenfreq",10));
	params->filefreq = int(inifile.getiniDouble("filefreq",10));
	params->thistfreq = int(inifile.getiniDouble("thistfreq",1));
	
	params->potparam1 = inifile.getiniDouble("potparam1",0.0);
	
	params->OutDir = inifile.getiniString("OutDir","output/");
	params->RunID = inifile.getiniString("RunID","run_01");
    checkdirexists(params->OutDir);
	
	params->flag = 0;
	
} // END GetParams()

void checkdirexists(string dir){
    
    using namespace boost::filesystem;
    
    if (!exists(dir + "/")) {
        cout << endl;
        cout << " --> Creating output directory" << endl;
        cout << endl;
        
        create_directory(dir);
    }
    
} // END checkdirexists()
