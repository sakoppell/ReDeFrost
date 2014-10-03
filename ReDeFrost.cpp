//
//  ReDeFrost.cpp
//  
//
//  Created by Stewart Koppell on 5/16/14.
//
//

#include <algorithm>
#include <complex>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

static const double PI=3.14159265359;
static const double x_max=10;
static const double sim_grid_size=128;
static const double dx=x_max/sim_grid_size;
static const int Np=sim_grid_size+4;//71;//primes: 64+4->71, 128+4->137, 256+4->263
static const int N=sim_grid_size+4;//sim_grid_size-4;
static const int Nscalars=2;
static const double lambda=.5;
static const double dt=dx*lambda;
static const int NNN=Np*Np*Np;


class field_class{
    public:
        double *field, *G, *Ginv, *M, *Lh, *Ah;
};


class cons_class{
    public:
        int wait,length;
        double *Error, *Pressure, *DxMomentum, *DtDensity, *Density, *Phi_cons;
};

class out_class{
    public:
        int out_bock,fft_samples, wait, length;
        double *meanfield, *stdevfield, *fftwfield;
    
};

inline int cp(int x, int y, int z){//field #, conjugate, time, x,y,z
    //if(x>N || x<0 || y>N || y<0 || z>N || z<0) cout<<"cp overflow "<<x<<" "<<y<<" "<<z<<endl;
    return z+Np*(y+Np*x);
}

inline int cg(int mu, int nu, int x,int y, int z){
    return z+N*(y+N*(x+N*(mu+4*nu)));
}

void wrap(double field[]){
    int block=0;
    for(int H=0;H<Nscalars;H++){for(int C=0;C<2;C++){
        block=(H*2+C)*NNN;
        for(int I=0;I<N;I++){for(int J=0;J<N;J++){
            field[cp(0,I,J)+block]=field[cp(N-4,I,J)+block];
            field[cp(1,I,J)+block]=field[cp(N-3,I,J)+block];
            field[cp(N-2,I,J)+block]=field[cp(2,I,J)+block];
            field[cp(N-1,I,J)+block]=field[cp(3,I,J)+block];
        
            field[cp(I,0,J)+block]=field[cp(I,N-4,J)+block];
            field[cp(I,1,J)+block]=field[cp(I,N-3,J)+block];
            field[cp(I,N-2,J)+block]=field[cp(I,2,J)+block];
            field[cp(I,N-1,J)+block]=field[cp(I,3,J)+block];
            
            field[cp(I,J,0)+block]=field[cp(I,J,N-4)+block];
            field[cp(I,J,1)+block]=field[cp(I,J,N-3)+block];
            field[cp(I,J,N-2)+block]=field[cp(I,J,2)+block];
            field[cp(I,J,N-1)+block]=field[cp(I,J,3)+block];
    }}}}
}

double DetG(field_class Phi, int x, int y, int z){
    double Detsum=0;
    for(int I=0;I<4;I++){for(int J=0;J<3;J++){for(int K=0;K<2;K++){int L=(K+1)%2;
                Detsum+=pow(-1,(I+J+K)%2)*Phi.G[cg(I,0,x,y,z)]*Phi.G[cg(J+(J>=I),1,x,y,z)]*Phi.G[cg(K+(K>=J)+((K+K>=J)>=I),2,x,y,z)]*Phi.G[cg(L+(L>=J)+((L+L>=J)>=I),3,x,y,z)];
    }}}
    return sqrt(-Detsum);
}


double Dx(field_class Phi, int fi, int c, int x, int y, int z, int index){//fi=field index add "rank" argument
    if(index==0 && c==0){return Phi.field[NNN*(2*fi+1)+cp(x,y,z)];}
    if(index==0 && c==1){cout<<"problem in Dx!";return 0;}
    else{
    int coo[3];coo[0]=0;coo[1]=0;coo[2]=0;
    coo[index-1]=1;
    int block=NNN*(2*fi+c);
    return (-Phi.field[block+cp(x+2*coo[0],y+2*coo[1],z+2*coo[2])]+8*Phi.field[block+cp(x+1*coo[0],y+1*coo[1],z+1*coo[2])]-8*Phi.field[block+cp(x-1*coo[0],y-1*coo[1],z-1*coo[2])]+Phi.field[block+cp(x-2*coo[0],y-2*coo[1],z-2*coo[2])])/(12*dx);
    }
}

double Dx(double field[], int fi, int c, int x, int y, int z, int index){//fi=field index add "rank" argument
    int coo[3];coo[0]=0;coo[1]=0;coo[2]=0;
    coo[index-1]=1;
    int block=NNN*(2*fi+c);
    return (-field[block+cp(x+2*coo[0],y+2*coo[1],z+2*coo[2])]+8*field[block+cp(x+1*coo[0],y+1*coo[1],z+1*coo[2])]-8*field[block+cp(x-1*coo[0],y-1*coo[1],z-1*coo[2])]+field[block+cp(x-2*coo[0],y-2*coo[1],z-2*coo[2])])/(12*dx);
}

double Dxx(double field[], int fi, int c, int x, int y, int z, int index){
    int coo[3];coo[0]=0;coo[1]=0;coo[2]=0;
    coo[index-1]=1;
    int block=NNN*(2*fi+c);
    return (-field[block+cp(x+2*coo[0],y+2*coo[1],z+2*coo[2])]+16*field[block+cp(x+1*coo[0],y+1*coo[1],z+1*coo[2])]-30*field[block+cp(x,y,z)]+16*field[block+cp(x-1*coo[0],y-1*coo[1],z-1*coo[2])]-field[block+cp(x-2*coo[0],y-2*coo[1],z-2*coo[2])])/(12*dx*dx);
}

double Gradient(double field[], int fi, int c, int x, int y, int z){
    int block=NNN*(2*fi+c);
    return (-field[block+cp(x+2,y,z)]-field[block+cp(x,y+2,z)]-field[block+cp(x,y,z+2)]
            +16*(field[block+cp(x+1,y,z)]+field[block+cp(x,y+1,z)]+field[block+cp(x,y,z+1)]+field[block+cp(x-1,y,z)]+field[block+cp(x,y-1,z)]+field[block+cp(x,y,z-1)])
            -90*(field[block+cp(x,y,z)])
            -field[block+cp(x-2,y,z)]-field[block+cp(x,y-2,z)]-field[block+cp(x,y,z-2)]
            )/(12*dx*dx);
}

double Pot(field_class Phi, int x, int y, int z){
    return .5*pow(Phi.field[cp(x,y,z)],2)+.5*10000*pow(Phi.field[cp(x,y,z)],2)*pow(Phi.field[cp(x,y,z)+2*NNN],2);
}

double DPot(field_class Phi, int x, int y, int z, int wrt){//derivative is with respect to field wrt
    if (wrt==0) return Phi.field[cp(x,y,z)]+10000*Phi.field[cp(x,y,z)]*pow(Phi.field[cp(x,y,z)+2*NNN],2);
    else return 10000*Phi.field[cp(x,y,z)+NNN*2]*pow(Phi.field[cp(x,y,z)],2);
}


double Christoffel(field_class Phi, int I, int J, int K){//GAMMA^I_J_K
    if(I==0 & J==K) return pow(Phi.Ah[0],2)/Phi.Lh[0];
    else if(J==0 & I==K) return 1.0/Phi.Lh[0];
    else if(K==0 & I==J) return 1.0/Phi.Lh[0];
    else return 0;
}


double _Tuv(field_class Phi, int u, int v, int x, int y, int z){
    double T=0;
    double Lag=Pot(Phi,x,y,z);
    for(int H=0;H<Nscalars;H++){
        for(int I=0;I<4;I++){
            for(int J=0;J<4;J++){
                Lag+=Dx(Phi,H,0,x,y,z,I)*Dx(Phi,H,0,x,y,z,J)*Phi.Ginv[cg(I,J,x,y,z)]/2;
            }
        }
    }
    for(int I=0;I<Nscalars;I++){
        T+=Dx(Phi,I,0,x,y,z,u)*Dx(Phi,I,0,x,y,z,v);
    }
    return T-Phi.G[cg(u,v,x,y,z)]*Lag;
}

void D_Phi(double K_Phi[], field_class Phi){
    for(int H=0;H<Nscalars;H++){
    #pragma omp parallel for
    for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){
        int block=2*NNN*H;
        K_Phi[cp(I,J,K)+block]=Phi.field[cp(I,J,K)+block+NNN];
        K_Phi[cp(I,J,K)+block+NNN]=-Phi.field[cp(I,J,K)+block+NNN]*3/Phi.Lh[0] -DPot(Phi,I,J,K,H)+Gradient(Phi.field,H,0,I,J,K)/pow(Phi.Ah[0],2);
    }}}}
    return;
}

void D_Lh(double K_Lh[], field_class Phi){
    
    double kin=0;
    double pot=0;
    
    #pragma omp parallel for reduction(+:pot,kin)
    for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++) {for(int K=2;K<N-2;K++){
        pot+=Pot(Phi,I,J,K);
        for(int H=0;H<Nscalars;H++){
        kin+=pow(Phi.field[cp(I,J,K)+(2*H+1)*NNN],2);
        }
    }}}
    
    K_Lh[0]=1+pow(Phi.Lh[0],2)*((kin-pot)/pow(N,3))/3;//0;//
    if(isnan(K_Lh[0]+abs(K_Lh[0]))) K_Lh[0]=1./0.;
    return;
}

void D_Ah(double K_Ah[], field_class Phi,int inflation_on){
    K_Ah[0]=Phi.Ah[0]/Phi.Lh[0];
    if(inflation_on==0) K_Ah=0;
    
    #pragma omp parallel for
    for(int I=0;I<N;I++){for(int J=0;J<N;J++) {for(int K=0;K<N;K++){
        Phi.G[cg(0,0,I,J,K)]=-1;
        Phi.Ginv[cg(0,0,I,J,K)]=-1;
        for(int L=1;L<4;L++){
            Phi.G[cg(L,L,I,J,K)]=pow(Phi.Ah[0],2);
            Phi.Ginv[cg(L,L,I,J,K)]=pow(Phi.Ah[0],-2);
    }}}}
    
    return;
}

void startgrid(field_class Phi, int N, int Nscalars, double H0, double phi0, double pi0, int inflation_on){
    //initial conditions
    for(int I=0;I<2*NNN*Nscalars;I++) Phi.field[I]=0;
    for(int I=0;I<NNN*16;I++) Phi.G[I]=0;
    
    Phi.M[0]=1;
    Phi.M[1]=0;
    
    fstream BD_file;
    BD_file.open("BD_file.txt", ios::in | ios::binary);
    double * Ttemp;
    Ttemp=new double[int(pow(N-4,3))*Nscalars];
    //BD_file.read(reinterpret_cast<char*>(Ttemp), sizeof(double)*gridsize);
    BD_file.close();
    
    for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++) {for(int K=2;K<N-2;K++){
        Phi.field[cp(I,J,K)]=phi0;//sin((I-2)*2*PI/(N-4));//exp(-64/pow(double(N),2)*(pow((I-N/2),2)+pow((J-N/2),2)));//1;//
        Phi.field[cp(I,J,K)+NNN]=pi0;//sin(J*2*pi/N)*sin(I*2*pi/N);pi0;//
        Phi.field[cp(I,J,K)+2*NNN]=0;
        Phi.field[cp(I,J,K)+3*NNN]=0;
    }}}
    wrap(Phi.field);

    for(int I=0;I<N;I++){for(int J=0;J<N;J++) {for(int K=0;K<N;K++){
            Phi.G[cg(0,0,I,J,K)]=-1;
            Phi.Ginv[cg(0,0,I,J,K)]=-1;
        for(int L=1;L<4;L++){
            Phi.G[cg(L,L,I,J,K)]=1;
            Phi.Ginv[cg(L,L,I,J,K)]=1;
    }}}}

    Phi.Ah[0]=1;
    
    if(inflation_on==0) Phi.Lh[0]=1.0/0.0;
    else Phi.Lh[0]=1./H0;


    for(int I=0;I<N;I++){for(int J=0;J<N;J++) {for(int K=0;K<N;K++){
        Phi.G[cg(0,0,I,J,K)]=-1;
        Phi.Ginv[cg(0,0,I,J,K)]=-1;
        for(int L=1;L<4;L++){
            Phi.G[cg(L,L,I,J,K)]=pow(Phi.Ah[0],2);
            Phi.Ginv[cg(L,L,I,J,K)]=pow(Phi.Ah[0],-2);
    }}}}
    
    delete [] Ttemp;
    return;
}

void advance(field_class Phi,int inflation_on){

    int gridsize=2*NNN*Nscalars;
    
    field_class Phi_1;
    
    Phi_1.field=new double[gridsize];
    Phi_1.G=new double[16*NNN];
    Phi_1.Ginv=new double[16*NNN];
    Phi_1.Lh=new double[1];
    Phi_1.Ah=new double[1];
    
    double *K1_Phi, *K2_Phi, *K3_Phi, *K4_Phi;
    
    K1_Phi=new double[gridsize];
    K2_Phi=new double[gridsize];
    K3_Phi=new double[gridsize];
    K4_Phi=new double[gridsize];
    
    double K1_Lh[1], K2_Lh[1], K3_Lh[1], K4_Lh[1];
    double K1_Ah[1], K2_Ah[1], K3_Ah[1], K4_Ah[1];
    
    D_Phi(K1_Phi,Phi);
    D_Lh(K1_Lh,Phi);
    D_Ah(K1_Ah,Phi,inflation_on);
    
    #pragma omp parallel for
    for(int I=0;I<gridsize;I++){
        Phi_1.field[I]=Phi.field[I]+dt/2*K1_Phi[I];
    }
    
    Phi_1.Lh[0]=Phi.Lh[0]+dt/2*K1_Lh[0];
    Phi_1.Ah[0]=Phi.Ah[0]+dt/2*K1_Ah[0];
    
    wrap(Phi_1.field);
    
    D_Phi(K2_Phi,Phi_1);
    D_Lh(K2_Lh,Phi_1);
    D_Ah(K2_Ah,Phi_1,inflation_on);
    
    #pragma omp parallel for
    for(int I=0;I<gridsize;I++){
        Phi_1.field[I]=Phi.field[I]+dt/2*K2_Phi[I];
    }
    
    Phi_1.Lh[0]=Phi.Lh[0]+dt/2*K2_Lh[0];
    Phi_1.Ah[0]=Phi.Ah[0]+dt/2*K2_Ah[0];
    
    wrap(Phi_1.field);
    
    D_Phi(K3_Phi,Phi_1);
    D_Lh(K3_Lh,Phi_1);
    D_Ah(K3_Ah,Phi_1,inflation_on);
    
    #pragma omp parallel for
    for(int I=0;I<gridsize;I++){
        Phi_1.field[I]=Phi.field[I]+dt*K3_Phi[I];
    }
    
    Phi_1.Lh[0]=Phi.Lh[0]+dt*K3_Lh[0];
    Phi_1.Ah[0]=Phi.Ah[0]+dt*K3_Ah[0];

    wrap(Phi_1.field);

    D_Phi(K4_Phi,Phi_1);
    D_Lh(K4_Lh,Phi_1);
    D_Ah(K4_Ah,Phi_1,inflation_on);
    
    #pragma omp parallel for
    for(int I=0;I<gridsize;I++){
        Phi.field[I]+=dt/6*(K1_Phi[I]+2*K2_Phi[I]+2*K3_Phi[I]+K4_Phi[I]);
    }

    Phi.Lh[0]+=dt/6*(K1_Lh[0]+2*K2_Lh[0]+2*K3_Lh[0]+K4_Lh[0]);
    Phi.Ah[0]+=dt/6*(K1_Ah[0]+2*K2_Ah[0]+2*K3_Ah[0]+K4_Ah[0]);
    
    wrap(Phi.field);
    
    delete [] K1_Phi;
    delete [] K2_Phi;
    delete [] K3_Phi;
    delete [] K4_Phi;
    delete [] Phi_1.field;
    delete [] Phi_1.G;
    delete [] Phi_1.Ginv;
    delete [] Phi_1.Ah;
    delete [] Phi_1.Lh;
}

/*
void CheckConservation(cons_class Cons, fstream& grid_file){
    //"Tuv" now stores:
    //(U,V)=(0,0) the integral of the S-E tensor scaled by sqrt(-g) (at times t%Eblock<5)
    //(U,V)=(0,i) the integral of the derivative wrt i of the S-E tensor scaled by sqrt(-g) (at times t/Eblock==2)
    double * DtT00;
    DtT00=new double[Cons.length];
    //take time derivative with central difference method (4rth order accuracy)
    for(int I=0;I<Cons.length;I++){
        DtT00[I]=(-Cons.Tuv[16*(5*I+4)]+8*Cons.Tuv[16*(5*I+3)]-8*Cons.Tuv[16*(5*I+1)]+Cons.Tuv[16*(5*I)])/(12*dt);
    }
    
    for(int I=0;I<Cons.length;I++) Cons.Phi_cons[I]=(-Cons.Phi_cons[I+4]+8*Cons.Phi_cons[I+3]-8*Cons.Phi_cons[I+1]+Cons.Phi_cons[I])/(12*dt);
    
    
    for(int I=0;I<Cons.length;I++) Cons.Momentum[I]=Cons.Tuv[1+16*(I*5+2)]+Cons.Tuv[2+16*(I*5+2)]+Cons.Tuv[3+16*(I*5+2)];
    
    double * cons_E2;
    cons_E2=new double[Cons.length];
    
    for(int I=0;I<Cons.length;I++) cons_E2[I]=pow((DtT00[I]+Cons.Momentum[I]+3*Cons.Tuv[16*(I*5+2)]/Cons.Lhist[I]+Cons.Pressure[I]*pow(Cons.Ahist[I],2)/Cons.Lhist[I]),2);
    
    for(int I=0;I<Cons.length;I++){
        grid_file<<DtT00[I]<<"   "<<Cons.Momentum[I]<<"   "<<3*Cons.Tuv[16*(I*5+2)]/Cons.Lhist[I]<<"  "<<Cons.Pressure[I]/pow(Cons.Ahist[I],2)/Cons.Lhist[I]<<"   "
                 <<(DtT00[I]+Cons.Momentum[I]+3*Cons.Tuv[16*(I*5+2)]/Cons.Lhist[I]+Cons.Pressure[I]*pow(Cons.Ahist[I],2)/Cons.Lhist[I])<<endl;
        grid_file<<sqrt(cons_E2[I])/NNN<<endl;
    }
    grid_file<<endl<<endl;

    
    for(int I=0;I<Cons.length;I++) cons_E2[0]+=cons_E2[I];

    grid_file<<"sqrt(E^2) per grid point "<<sqrt(cons_E2[0]/Cons.length)/pow(N,3)<<endl;
    Cons.Phi_cons[0]=pow(Cons.Phi_cons[0],2);

    for(int I=0;I<Cons.length;I++){Cons.Phi_cons[0]+=pow(Cons.Phi_cons[I],2);}
    
    grid_file<<"sqrt(phi^2) per grid point "<<sqrt(Cons.Phi_cons[0]/Cons.length)/pow(N,3)<<endl;
    
    grid_file<<"residual curvature "<<endl;
    for(int I=0;I<Cons.length;I++) grid_file<<Cons.Ahist[I]*Cons.Ahist[I]*(Cons.Tuv[16*(5*I+2)]/3-1./Cons.Lhist[I]/Cons.Lhist[I])<<endl;
 
    delete [] cons_E2;
    delete [] DtT00;
    return;
}
*/

void Ecollect(field_class Phi, cons_class Cons, int t, int Tmax){
    if(t%(Cons.wait+5)==1){
        double pressure=0;
        int E_t=t/(Cons.wait+5);
        #pragma omp parallel for reduction(+:pressure)
        for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){for(int L=1;L<4;L++){
                        pressure+=_Tuv(Phi,L,L,I,J,K);
        }}}}
        
        Cons.Pressure[E_t]=pressure/pow(N-4,3)/3*pow(Phi.Ah[0],2);
        
        double * Ttemp;
        Ttemp=new double[NNN];
        double dxmomentum=0;
        
        for(int U=1;U<4;U++){
        #pragma omp parallel for
        for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){
                    Ttemp[cp(I,J,K)]=_Tuv(Phi,U,0,I,J,K);
        }   }   }
        dxmomentum=0;
        #pragma omp parallel for reduction(+:dxmomentum)
        for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){
                    dxmomentum+=Dx(Ttemp,0,0,I,J,K,U);
        }   }   }
        }
        delete [] Ttemp;
        Cons.DxMomentum[E_t]=dxmomentum/pow(N-4,3)/3;
    }
    
    int E_t=t/(Cons.wait+5)*5+t%(Cons.wait+5);//Current E index
    
    double density=0;
    
    #pragma omp parallel for reduction(+:density)
    for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){
                density+=_Tuv(Phi,0,0,I,J,K);
    }   }   }
    
    Cons.Density[E_t]=density/pow(N-4,3);
  
    if(t%(Cons.wait+5)==4){
        int E_n=t/(Cons.wait+5);
        Cons.DtDensity[E_n]=(Cons.Density[E_t-4]-8*Cons.Density[E_t-3]+8*Cons.Density[E_t-1]-Cons.Density[E_t])/(12*dt);
        Cons.Error[E_n]=Cons.DtDensity[E_n]+3*Cons.DxMomentum[E_n]+3.0/Phi.Lh[0]*(pow(Phi.Ah[0],2)*Cons.Pressure[E_n]+Cons.Density[E_t-2]);

    //cout<<Cons.Density[E_n]<<endl;
    //cout<<Cons.Pressure[E_n]<<endl;
    }
    
    
    /*for(int U=0;U<4;U++){
    #pragma omp parallel for
    for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){
                Cons.Phi_cons[E_t]+=DetG(Phi,I,J,K)*Phi.Ginv[cg(0,U,I,J,K)]*Dx(Phi.field,0,0,I,J,K,U);
                if(U==0) Cons.Phi_cons[E_t+5*Cons.length]=DetG(Phi,I,J,K)*DPot(Phi,I,J,K,0);
    }   }   }
    }*/
    
    return;
}

void Ocollect(field_class Phi, out_class Fout, int T, fstream& grid_file){
    double temp=0;
    for(int H=0;H<Nscalars*2;H++){
        temp=0;
        #pragma omp parallel for reduction(+:mean)
        for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){ temp+=Phi.field[cp(I,J,K)+H*NNN];}}}
        Fout.meanfield[T/Fout.wait+H*Fout.length]=temp/pow(N-4,3);
    }
    for(int H=0;H<Nscalars*2;H++){
        temp=0;
        #pragma omp parallel for reduction(+:mean)2*H*
        for(int I=2;I<N-2;I++){for(int J=2;J<N-2;J++){for(int K=2;K<N-2;K++){ temp+=pow(Phi.field[cp(I,J,K)+H*NNN]-Fout.meanfield[T/Fout.wait+H*Fout.length],2);}}}
        Fout.stdevfield[T/Fout.wait+H*Fout.length]=sqrt(temp)/pow(N-4,3);
    }

    //for(int H=0;H<Nscalars;H++) grid_file<<"*** Field "<<H<<" - mean phi: "<<Fout.meanfield[T/O_wait+2*H*OL]<<", stdev phi: "<<Fout.stdevfield[T/O_wait+2*H*OL]<<" - mean pi: "<<Fout.meanfield[T/O_wait+OL+2*H*OL]<<", stdev pi: "<<Fout.stdevfield[T/O_wait+OL+2*H*OL]<<endl;
}

void collect(field_class Phi, out_class Fout, cons_class Cons, int t, int output_type, fstream& grid_file){

    if(output_type==1){
        if(t==0){
            grid_file.write(reinterpret_cast<const char*>(&N), sizeof(N));
            grid_file.write(reinterpret_cast<const char*>(&Nscalars), sizeof(Nscalars));
            grid_file.write(reinterpret_cast<const char*>(&dt), sizeof(dt));
        }
        grid_file.write(reinterpret_cast<const char*>(Phi.field), sizeof(double)*N*N*N*Nscalars);//reinterpret_cast<const char*>(&Phi[I])
    }
    if(output_type==2){
        grid_file<<t<<" "<<t*dt;
        //for(int I=2;I<N-2;I=I+(N-4)/1) grid_file<<"   "<<Phi.field[cp(I,2,2)]*pow(Phi.Ah[0],1.5);
        //for(int I=2;I<N-2;I=I+(N-4)/1) grid_file<<"   "<<Phi.field[cp(I,2,2)+NNN]*pow(Phi.Ah[0],1.5);
        for(int H=0;H<Nscalars;H++) grid_file<<"    "<<Fout.meanfield[t/Fout.wait+2*Fout.wait*H]*pow(Phi.Ah[0],1.5)<<"   "<<Fout.stdevfield[t/Fout.wait+2*Fout.wait*H]*pow(Phi.Ah[0],1.5)<<"  "<<Fout.meanfield[t/Fout.wait+2*Fout.wait*H+Fout.length]*pow(Phi.Ah[0],1.5)<<"   "<<Fout.stdevfield[t/Fout.wait+2*Fout.wait*H+Fout.length]*pow(Phi.Ah[0],1.5);
        grid_file<<"    "<<Phi.Ah[0]<<"    "<<Phi.Lh[0];
        
        if(t%(Cons.wait+5)==4) grid_file<<"  -  "<<Cons.DtDensity[t/(Cons.wait+5)]<<"  "<<Cons.DxMomentum[t/(Cons.wait+5)]<<"   "<<Cons.Pressure[t/(Cons.wait+5)]<<"    "<<Cons.Density[t/(Cons.wait+5)*5+2]<<" "<<Cons.Error[t/(Cons.wait+5)];
        grid_file<<endl;
    }

    return;
}

int main(){
    
    double t_max=10*dt;
    int datablock=1;
    int E_wait=1;
    int O_wait=1;
    int inflation_on=1;//1=on, 0=off
    int output_type=2;
    double H0=.50467;
    double phi0=1.009343;
    double pi0=-phi0/sqrt(2)*inflation_on;
    
    
    int Tmax=int(t_max/dt);
    int gridsize=2*Np*Np*Np*Nscalars;

    field_class PHI;
    PHI.field=new double[gridsize];
    PHI.G=new double[16*NNN];
    PHI.Ginv=new double[16*NNN];
    PHI.M=new double[Nscalars];
    PHI.Ah=new double[1];
    PHI.Lh=new double[1];

    cons_class CONS;

    int EL=(Tmax+1)/(E_wait+5);
    CONS.length=EL;
    CONS.wait=E_wait;
    CONS.Error=new double[EL];
    CONS.Pressure=new double[EL];
    CONS.DxMomentum=new double[EL];
    CONS.DtDensity=new double[EL];
    CONS.Density=new double[EL*5];
    CONS.Phi_cons=new double[EL*5];
    
    out_class FOUT;
    
    int OL=(Tmax)/(O_wait)+1;
    FOUT.wait=O_wait;
    FOUT.length=OL;
    FOUT.meanfield=new double[OL*Nscalars*2];
    FOUT.stdevfield=new double[OL*Nscalars*2];
    FOUT.fftwfield=new double[OL*20*Nscalars*2];
    
    
    fstream grid_file;
    grid_file.open("grid_file.txt", ios::out | ios::binary);
    
    startgrid(PHI,N,Nscalars,H0,phi0,pi0,inflation_on);
    
    time_t start_timer;
    time_t end_timer;
    time(&start_timer);
    
    Ocollect(PHI, FOUT, 0, grid_file);
    Ecollect(PHI, CONS, 0, Tmax);
    collect(PHI, FOUT, CONS, 0, output_type, grid_file);

    cout<<"0%  ---*"<<endl;
    for(int T=1;T<Tmax+1;T++){
        cout<<T*1.0/Tmax*100<<"%";
        advance(PHI,inflation_on);
        if (T%(E_wait+5)<5 && EL>0){ Ecollect(PHI, CONS, T, Tmax); cout<<"  ---";}
        if ((T)%(O_wait)==0){ Ocollect(PHI, FOUT, T, grid_file); cout<<" *";}
        if (T%O_wait==0) collect(PHI, FOUT, CONS, T, output_type, grid_file);
        cout<<endl;
    }
    time(&end_timer);
    cout<<"Main Loop Runtime: "<<difftime(end_timer,start_timer)<<endl;
    
    grid_file.close();
    return 1;
}