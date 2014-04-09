

// evolve.cpp

int Vi(int nn,int i, int j, int k, int *Vd){
    return nn*Vd[1]*Vd[2]*Vd[3]+i*Vd[2]*Vd[3]+j*Vd[3]+k;
}
int Si(int i, int j, int k, int *Sd){
    return  i*Sd[1]*Sd[2]+j*Sd[2]+k;
}


void evolve(){
    
    grid.h2=2.0*grid.h;
    grid.hh=pow(grid.h,2.0);
    grid.htht=pow(grid.ht,2.0);
    
    for(int t=1; t<tmax; t++){// loop over time-steps
        
        tt=1; tp=0;
        if(t%2==0){tt=0;tp=1;}
        
        a = GetScaleFactor(t);
        SolvePoisson();
        RunOverGrid();
        
    } // end t-loop

} // end evolve()

struct gridproto{
    
    double h;
    double h2;
    double hh;
    double ht;
    double htht;
    
    int imax;
    int jmax;
    int kmax;
    
    int i;
    int j;
    int k;
    int ip;
    int jp;
    int kp;
    int im;
    int jm;
    int km;
    
};

void RunOverGrid(){
    
    for(int i=0;i<grid.imax; i++){
        // Get i+ & i- (takes into account boundary condition)
        grid.i = i;
        grid.ip=GetCoord_p(i,imax);
        grid.im=GetCoord_m(i,imax);
        for(int j=0;j<grid.jmax;j++){
            // Get j+ & j- (takes into account boundary condition)
            grid.j = i;
            grid.jp=GetCoord_p(j,jmax);
            grid.jm=GetCoord_m(j,jmax);
            for(int k=0;k<grid.kmax;k++){
                // Get k+ & k- (takes into account boundary condition)
                grid.k = i;
                grid.kp=GetCoord_p(k,imax);
                grid.km=GetCoord_m(k,imax);
                
                // Call routine to do things at this grid location
                DoOnGrid(tt,&grid,field,)
                
            }// end k-loop
        }// end j-loop
    } // end i-loop
    
} // END RunOverGrid

void SolvePoisson(){
    
    double GSaccuracy,GSerror;
    int n;

    // Solve Poisson equation via Gauss-Siedel or SoR

    cout << "starting Gauss-Siedel solution of Poisson equation" << endl;
    
    GSstep=0;
    while(true){
        GSerror = 0.0;
        n = 0;
        nn = 1; np = 0;
        if(GSstep%2==0){nn=0;np=1;}
        for(int i=iGSmin;i<=iGSmax;i++){
            ip = i+1;
            im = i-1;
            for(int j=jGSmin;j<=jGSmax;j++){
                jp = j+1;
                jm = j-1;
                for(int k=kGSmin;k<=kGSmax;k++){
                    kp = k+1;
                    km = k-1;
                    V[Vi(np,i,j,k)]=(V[Vi(nn,ip,j,k)]+V[Vi(nn,i,jp,k)]+V[Vi(nn,i,j,kp)]+V[Vi(np,im,j,k)]+V[Vi(np,i,jm,k)]+V[Vi(np,i,j,km)]-h2*S[Si(i,j,k)])/6.0;
                    // Compute GS error
                    if(V[Vi(np,i,j,k)]!=V[Vi(nn,i,j,k)]){
                        GSerror+=abs(1.0-V[Vi(np,i,j,k)]/V[Vi(nn,i,j,k)]);
                        n++;
                    }
                    
                } // END k-loop
            } // END j-loop
        } // END i-loop
        
        
        if(n!=0) GSerror=GSerror/n;
        
        // Stop Gauss-Sidel only when error less than accuracy
        if(GSerror < GSaccuracy)
            break;
        
    } // END Gauss-Siedel solve of Poisson equation
    
} // END SolvePoisson




void DoOnGrid(){

    // Do these things at this grid site
    
    // Compute derivatives of the field
    GetFieldDerivs();
    
    // Construct equation of motion
    // this can be inside the previous c-loop
    // if no derivatives of field c' are in
    // Eom of field c. Keep them separate for now,
    // for generality
    
    for(int c=0; c<cmax;c++){
        // This is an example for the flat-space KG equation
        lap = derivs[ derivsi(c,3) ] + derivs[ derivsi(c,4) ] + derivs[ derivsi(c,5) ];
        eom[c] = lap - dpot[c];
        
    }// end c-loop
    
    // Update field value
    UpdateField();
    
    // Get the source for the Poisson equation
    dum=0.0;
    for(int c=0;c<cmax;c++){
        dum=dum+pow(field[tt][c][i][j][k],2.0);
    }
    
    PoissSource[i][j][k]=4.0*pi*G*rho0*(dum-1.0)/a;
    
} // END DoOnGrid

void UpdateField(){
    
    for(int c=0;c<cmax;c++){
        field[tp][c][i][j][k]=eom[c]*htht+2.0*field[tt][c][i][j][k]-field[tp][c][i][j][k];
    }
    
} // END UpdateField

void GetFieldDerivs(){
    
    if(derivsaccuracy==2){
        // Second order accurate derivatives
        
        for(int c=0;c<cmax; c++){
            // Get field forwards & backwards from current location
            fijk=field[tt][c][i][j][k];
            fip=field[tt][c][ip][j][k];
            fim=field[tt][c][im][j][k];
            fjp=field[tt][c][i][jp][k];
            fjm=field[tt][c][i][jm][k];
            fkp=field[tt][c][i][j][kp];
            fkm=field[tt][c][i][j][km];
            
            // First derivatives
            derivs[ derivsi(c,0) ] = (fip-fim)/h2;
            derivs[ derivsi(c,1) ] = (fjp-fjm)/h2;
            derivs[ derivsi(c,2) ] = (fkp-fkm)/h2;
            // Second derivatives
            derivs[ derivsi(c,3) ] = (fip+fim-2.0*fijk)/hh;
            derivs[ derivsi(c,4) ] = (fjp+fjm-2.0*fijk)/hh;
            derivs[ derivsi(c,5) ] = (fkp+fkm-2.0*fijk)/hh;

        } // end c-loop
        
    }
    if(derivsaccuracy==4){
        // Fourth order accurate derivatives
    }

} // END derivs

int derivsi(int c, int ID){
 
    return c*cmax+ID;
    
}

int GetCoord_p(int i,int imax){
    
    ip=i+1
    if(i==imax-1){ip=0;}
    return ip;
    
} // END GetCoord_p

int GetCoord_m(int i,int imax){
    
    im=i-1;
    if(i==0){im=imax-1;}
    return im;
    
} // END GetCoord_m

double GetScaleFactor(int t){
 
    power = 2.0/3.0;
    return a = pow(t*ht,power);
    
} // end getscalefactor()
