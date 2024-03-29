
// setup.cpp

// This contains routines to read in the params.ini file
//  and populate items in the structs

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
	
} // END SetupField()

void SetupLaplacianStencil(struct DATA *params, struct LAPLACIANSTENCIL *stencil){

	int whichstencil = params->lapstencilchoice;
	stencil->lapstencilchoice = whichstencil;
	
	if( whichstencil == 1 ){
		stencil->c0 = 6.0;
		stencil->c1 = 1.0;
		stencil->c2 = 0.0;
		stencil->c3 = 0.0;
	}
	
	if( whichstencil == 2 ){
		stencil->c0 = 4.0;
		stencil->c1 = 1.0/3.0;
		stencil->c2 = 1.0/6.0;
		stencil->c3 = 0.0;
	}
	
	if( whichstencil == 3 ){
		stencil->c0 = 14.0/3.0;
		stencil->c1 = 2.0/3.0;
		stencil->c2 = 0.0;
		stencil->c3 = 1.0/12.0;
	}
	
	if( whichstencil == 4 ){
		stencil->c0 = 64.0/15.0;
		stencil->c1 = 7.0/15.0;
		stencil->c2 = 1.0/10.0;
		stencil->c3 = 1.0/30.0;
	}

} // END SetupLaplacianStencil()

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
    params->fx_i = int(inifile.getiniDouble("fix_i",0));
    params->fx_j = int(inifile.getiniDouble("fix_j",0));
    params->fx_k = int(inifile.getiniDouble("fix_k",0));    
	params->pottype = int(inifile.getiniDouble("pottype",1));
	params->inittype = int(inifile.getiniDouble("inittype",1));
	params->evoltype = int(inifile.getiniDouble("evoltype",0));
	params->eomtype = int(inifile.getiniDouble("eomtype",0));
	params->lapstencilchoice = int(inifile.getiniDouble("lapstencilchoice",1));
	params->screenfreq = int(inifile.getiniDouble("screenfreq",10));
	params->filefreq = int(inifile.getiniDouble("filefreq",10));
	params->thistfreq = int(inifile.getiniDouble("thistfreq",1));
	params->potparam1 = inifile.getiniDouble("potparam1",0.0);
	params->OutDir = inifile.getiniString("OutDir","output/");
	params->RunID = inifile.getiniString("RunID","run_01");
	
    checkdirexists(params->OutDir);
	
	params->flag = 0;
	
	
	
} // END GetParams()

void CheckParams(struct DATA *params){

	// This checks the parameters for sanity reasons
	// If anything has been setup pathologically, code will not run
	
	if( params->imax < 0 ) params->flag = 1;
	if( params->jmax < 0 ) params->flag = 1;
	if( params->kmax < 0 ) params->flag = 1;
	
	if( params->fx_i > params->imax ) params->flag = 1;
	if( params->fx_j > params->jmax ) params->flag = 1;
	if( params->fx_k > params->kmax ) params->flag = 1;
	
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