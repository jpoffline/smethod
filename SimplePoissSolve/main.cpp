
#include "main.h"

const int imin = 0;
const int jmin = 0;
const int imax = 50;
const int jmax = 50;
const int nparams_iterate = 6;
double V[2][imax+1][jmax+1];
double S[imax+1][jmax+1];
double *params_iterate = new double[nparams_iterate];

string OutDir;
double h;
double accuracy;
int methodID;
double param1;

// Fixed x & y values to output at

double xf;
double yf;


int main(int argc, char* argv[]){
 
    
    // Get the parameter file:
    IniReader inifile;
    // If the user specified a param file at runtime
    //  then use that, else, use the default, params.ini
    //   (User specifies via ./EXE my_param.par)
    
    if (argc > 1)
        inifile.read(argv[1]);
    else
        inifile.read("params.ini");
    
    methodID = int(inifile.getiniDouble("methodID",1));
    h = inifile.getiniDouble("h",0.1);
    param1 = inifile.getiniDouble("param1",0.0);
    accuracy = inifile.getiniDouble("accuracy",0.0001);
    OutDir = inifile.getiniString("OutDir","output/");
    xf = inifile.getiniDouble("xf",0.0);
    yf = inifile.getiniDouble("yf",0.0);
    
    checkdirexists(OutDir);
    int test_number = 1;
    
    params_iterate[0] = imin+1;
    params_iterate[1] = imax-1;
    params_iterate[2] = jmin+1;
    params_iterate[3] = jmax-1;
    params_iterate[4] = h*h;
    params_iterate[5] = accuracy;
    
    /*
    params_init[0]=test_number;
    params_init[1]=param1;
    */
    // (1) Initialise the source & potential
    initialize_source_pot(test_number,h);
    
    // (2) Iterate using chosen method (GS or SoR)
    iterate(methodID);
    
    // (3) Output the final configuration
    output_finalconfig(methodID);
    
    cout << endl;
    delete params_iterate;
    
}// end main()


void initialize_source_pot(int test_number, double h){
    // Setup source & potential
    
    
    double softening,charge,dist,x,y;

    if(test_number==0){
        cout << "test = 0" << endl;
        softening=0.1;
        charge=10.0;
    }
    if(test_number==1){
        cout << "ic test = 1 " << endl;
        if(param1<1) cout << " > Source and potential initially exact" << endl;
        if(param1>2) cout << " > Source and potential not exact" << endl;
    }
    
    ofstream icout,ic_x,ic_y;
    icout.open(OutDir+"ic.dat");
    ic_x.open(OutDir+"ic_x.dat");
    ic_y.open(OutDir+"ic_y.dat");

    if(test_number==0){
        for(int i=0;i<=imax;i++){
            for(int j=0;j<=jmax;j++){
                
                x = (i-imax*0.5)*h;
                y = (j-jmax*0.5)*h;

                x = x+softening;
                y = y+softening;
                
                dist = sqrt(x*x+y*y);
                S[i][j] = charge/dist;
                V[0][i][j]=1.0;
                V[1][i][j]=V[0][i][j];
                icout << x << " " << y << " " << V[0][i][j] << " " << S[i][j] << endl;
                if(x==xf){
                    ic_x << x << " " << y << " " << V[0][i][j] << " " << S[i][j] << endl;
                }
                if(y==yf){
                    ic_y << x << " " << y << " " << V[0][i][j] << " " << S[i][j] <<  endl;
                }
                
            } // end j-loop
        } // end i-loop
    } // end if(test_number==0)
    
    if(test_number==1){

        
        for(int i=0;i<=imax;i++){
            for(int j=0;j<=jmax;j++){
                
                x = (i-imax*0.5)*h;
                y = (j-jmax*0.5)*h;
                S[i][j] = 2.0*pow(x,3.0)+6.0*x*y*(y-1.0);
                if(param1<1) V[0][i][j] = pow(x,3.0)*y*(1.0-y); // exact solution
                if(param1>2) V[0][i][j] = pow(x,3.0)*x*(1.0-y); // not exact

                V[1][i][j] = V[0][i][j];

                icout << x << " " << y << " " << V[0][i][j] << " " << S[i][j] << endl;
                if(x==xf){
                    ic_x << x << " " << y << " " << V[0][i][j] << " " << S[i][j] << endl;
                }
                if(y==yf){
                    ic_y << x << " " << y << " " << V[0][i][j] << " " << S[i][j] <<  endl;
                }
                
            } // end j-loop
            icout << endl;
        } // end i-loop
       
        
    } // end if(test_number==1)
    
    icout.close();
    ic_x.close();
    ic_y.close();
    
} // end initialize_source_pot()

void output_finalconfig(int methodID){
    
    double x,y;
    ofstream finalout,xfix,yfix;
    
    if(methodID==1){
        finalout.open(OutDir+"finalconfig_GS.dat");
        xfix.open(OutDir+"finalconfig_GS_xfix.dat");
        yfix.open(OutDir+"finalconfig_GS_yfix.dat");
    }
    if(methodID==2){
        finalout.open(OutDir+"finalconfig_SoR.dat");
        xfix.open(OutDir+"finalconfig_SoR_xfix.dat");
        yfix.open(OutDir+"finalconfig_SoR_yfix.dat");
    }
    
    cout << "xfix = " << xf << endl;
    cout << "yfix = " << yf << endl;
    
    for(int i=0;i<=imax;i++){
        for(int j=0;j<=jmax;j++){
            x = (i-imax*0.5)*h;
            y = (j-jmax*0.5)*h;
            finalout << x << " " << y << " " << V[0][i][j] << " " << S[i][j] << endl;

            if(x==xf){
                xfix << x << " " << y << " " << V[0][i][j] << " " << S[i][j] << endl;
            }
            if(y==yf){
                yfix << x << " " << y << " " << V[0][i][j] << " " << S[i][j] <<  endl;
            }
        } // end j-loop
        finalout << endl;
    } // end i-loop
    finalout.close();
    xfix.close();
    yfix.close();
    
} // end output_finalconfig()

void iterate(int methodID){
    
    // This keeps updating the potential,
    //  using GS or SoR.
    
    // Terminates when the error is less than
    //  desired accuracy.
    
    int step=0;
    double accuracy = params_iterate[5];
    double error;
    
    if(methodID == 1){
        cout << "Gauss-Siedel" << endl;
    }
    if(methodID == 2){
        cout << "SoR" << endl;
    }
    
    ofstream outerror;
    
    if(methodID==1) outerror.open(OutDir+"err_GS.dat");
    if(methodID==2) outerror.open(OutDir+"err_SoR.dat");
    
    while(true){
        
        if(methodID == 1){
            gauss_seidel(step);
        }
        if(methodID == 2){
            successive_over_relaxation(step);
        }
        step++;
        
        error = compute_error(step);
        if(error < accuracy)
            break;
        
        outerror << step << " " << error << endl;
        
    }
    outerror.close();
    
    cout << endl;
    cout << "Error less than accuracy" << endl;
    cout << "error = " << error <<  endl;
    cout << "accuracy = " << accuracy << endl;
    cout << "Number of steps = " << step << endl;
    
    
} // end iterate

double compute_error(int step){
    
    int is_min = int(params_iterate[0]);
    int is_max = int(params_iterate[1]);
    int js_min = int(params_iterate[2]);
    int js_max = int(params_iterate[3]);
    
    double error = 0;
    int n = 0;
    
    int nn = 1;
    int np = 0;
    if(step%2==0){nn=0; np=1;}
    
    double newpot,oldpot;
    
    for(int i=is_min;i<=is_max;i++){
        for(int j=js_min;j<=js_max;j++){
            newpot=V[np][i][j];
            
            if(newpot!=0){
                oldpot=V[nn][i][j];
                if(newpot!=oldpot){
                    error+=abs(1.0-newpot/oldpot);
                    ++n;
                } // end if()
            }// end if()
            
        } // end j-loop
    }// end i-loop
    
    if(n!=0) error = error/n;
    
    return error;
    
} // end compute_error()

void gauss_seidel(int step){
    
    double h2 = params_iterate[4];
    
    int is_min = int(params_iterate[0]);
    int is_max = int(params_iterate[1]);
    int js_min = int(params_iterate[2]);
    int js_max = int(params_iterate[3]);
    
    int nn = 1;
    int np = 0;
    if(step%2==0){nn=0; np=1;}
    
    for(int i=is_min;i<=is_max;i++){
        for(int j=js_min;j<=js_max;j++){
            
            V[np][i][j]=0.25*(V[nn][i+1][j]+V[np][i-1][j]+V[nn][i][j+1]+V[np][i][j-1]+h2*S[i][j]);
            
        } // end j-loop
    } // end i-loop
    
} // end gauss_seidel

void successive_over_relaxation(int step){
    
    double h2=params_iterate[4];
    
    int is_min = int(params_iterate[0]);
    int is_max = int(params_iterate[1]);
    int js_min = int(params_iterate[2]);
    int js_max = int(params_iterate[3]);
    
    // SoR parameter
    double s = 2.0/(1.0+pi/imax);
    
    int nn=1;
    int np=0;
    if(step%2==0){nn=0; np=1;}
    
    for(int i=is_min;i<=is_max;i++){
        for(int j=js_min;j<=js_max;j++){
            
            V[np][i][j]=(1.0-s)*V[nn][i][j]+0.25*s*(V[nn][i+1][j]+V[np][i-1][j]+V[nn][i][j+1]+V[np][i][j-1]+h2*S[i][j]);
            
        } // end j-loop
    } // end i-loop
    
} // end successive_over_relaxation

void checkdirexists(string dir){
    
    using namespace boost::filesystem;
    
    if (!exists(dir + "/")) {
        cout << endl;
        cout << " --> Creating output directory" << endl;
        cout << endl;
        
        create_directory(dir);
    }
    
} // END checkdirexists


