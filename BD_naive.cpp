//
//  BD_naive.cpp
//  BD_naive
//
//  Created by Stewart Koppell on 8/11/14.
//  Copyright (c) 2014 Stewart Koppell. All rights reserved.
//

#include "BD_naive.h"
#include "ReDeFrost2.h"

int N=int(x_max/dx);
int GS=N*N*N;
double dk=2*PI/x_max;

void grand(double GRAND[])//gausian random noise
{
    SEED++;

    srand(int(time(NULL))*SEED+SEED);
    int RANDSTRING=rand();
    URAND=(RANDSTRING%10000)/double(10000);
    VRAND=(((RANDSTRING-(RANDSTRING%10000+1))/10000)%10000+1)/double(10000);
    if(URAND==0){
        //cout<<"-------------- URAND=0 ---------------"<<endl;
        URAND=0.1426837;
    }
    if(VRAND==0){
        //cout<<"-------------- VRAND=0 ---------------"<<endl;
        VRAND=0.8932567;
    }
    if(rand_on==1){
        GRAND[0]=sqrt(-2.0*log(URAND))*cos(2*PI*VRAND);
        GRAND[1]=sqrt(-2.0*log(URAND))*sin(2*PI*VRAND);
    }else{
        GRAND[0]=cos(PI/4);
        GRAND[1]=sin(PI/4);
    }
    if (URAND==0 || VRAND==0) cout<<URAND<<" "<<VRAND<<endl;
    return;
}


int main(){
    Meff=new double [Nscalars];
    for(int H=0;H<Nscalars;H++) Meff[H]=9/4.*H0+M[H];
    phi_x=new fftw_complex[GS];
    phi_k=new fftw_complex[GS];
    pi_x=new fftw_complex[GS];
    pi_k=new fftw_complex[GS];
    BD_out=new double[GS*Nscalars*2];
    phi_out=new double[GS];
    symtrack=new int[GS];
    fftw_plan plan_phi;
    fftw_plan plan_pi;
    plan_phi=fftw_plan_dft_3d(N, N, N, phi_k, phi_x,FFTW_BACKWARD,FFTW_MEASURE);
    plan_pi=fftw_plan_dft_3d(N,N,N,pi_k,pi_x,FFTW_BACKWARD,FFTW_MEASURE);
    
    double kk=0;
    double rescale=1.0/sqrt(2*PI*dk*dk*dk)/sqrt(2);
    int kx=0;
    int ky=0;
    int kz=0;
    int coo=0;
    //int Acoo=0;
    for(int H=0;H<Nscalars;H++){
        for(int Kx=0;Kx<N;Kx++){for(int Ky=0;Ky<N;Ky++){for(int Kz=0;Kz<N;Kz++){
            kx=Kx;ky=Ky;kz=Kz;
            if(Kx>double(N)/2) kx=N-Kx;if(Ky>double(N)/2) ky=N-Ky;if(Kz>double(N)/2) kz=N-Kz;
            kk=(pow(kx,2)+pow(ky,2)+pow(kz,2))*dk*dk;
            coo=Kx+N*(Ky+N*Kz);
            phi_k[coo][0]=1./sqrt(sqrt(kk+Meff[H]*Meff[H]))*rescale;
            pi_k[coo][0]=sqrt(sqrt(kk+Meff[H]*Meff[H]))*rescale;
        }}}
        
        //enforce sym
        for(int I=0;I<GS;I++) symtrack[I]=1;
        double Grand[2];
        for(int Kx=0;Kx<N;Kx++){for(int Ky=0;Ky<N;Ky++){for(int Kz=0;Kz<N;Kz++){
            kx=(N-Kx)%N; ky=(N-Ky)%N; kz=(N-Kz)%N;
            if(kx==Kx && ky==Ky && kz==Kz && symtrack[Kx+N*(Ky+N*Kz)]==1){
                symtrack[kx+N*(ky+N*kz)]=0;
                grand(Grand);
                phi_k[Kx+N*(Ky+N*Kz)][0]=phi_k[Kx+N*(Ky+N*Kz)][0]*sqrt(Grand[0]*Grand[0]+Grand[1]*Grand[1]);
                phi_k[Kx+N*(Ky+N*Kz)][1]=0;
                grand(Grand);
                pi_k[Kx+N*(Ky+N*Kz)][0]=phi_k[Kx+N*(Ky+N*Kz)][0]*sqrt(Grand[0]*Grand[0]+Grand[1]*Grand[1]);
                pi_k[Kx+N*(Ky+N*Kz)][1]=0;
            }
            else if(symtrack[Kx+N*(Ky+N*Kz)]==1){
                symtrack[kx+N*(ky+N*kz)]=0;
                grand(Grand);
                phi_k[kx+N*(ky+N*kz)][1]=phi_k[kx+N*(ky+N*kz)][0]*Grand[1];
                phi_k[kx+N*(ky+N*kz)][0]=phi_k[kx+N*(ky+N*kz)][0]*Grand[0];
                grand(Grand);
                pi_k[kx+N*(ky+N*kz)][1]=phi_k[kx+N*(ky+N*kz)][0]*Grand[1];
                pi_k[kx+N*(ky+N*kz)][0]=phi_k[kx+N*(ky+N*kz)][0]*Grand[0];
            }
            else{
                phi_k[kx+N*(ky+N*kz)][0]=phi_k[Kx+N*(Ky+N*Kz)][0];
                phi_k[kx+N*(ky+N*kz)][1]=-phi_k[Kx+N*(Ky+N*Kz)][1];
                pi_k[kx+N*(ky+N*kz)][0]=phi_k[Kx+N*(Ky+N*Kz)][0];
                pi_k[kx+N*(ky+N*kz)][1]=-phi_k[Kx+N*(Ky+N*Kz)][1];
            }
        }}}

        //for(int I=0;I<GS;I++) phi_out[I]=pow(phi_k[I][0],2)+pow(phi_k[I][1],2);
    
        fftw_execute(plan_phi);
        fftw_execute(plan_pi);
        
        //for(int I=0;I<GS;I++) cout<<phi_x[I][0]<<"   "<<phi_x[I][1]<<endl;
        
        /*for(int I=0;I<GS;I++) phi_out[I]=phi_x[I][0];
        
        diag_file.open("diag_file.txt", ios::out | ios::binary | ios::trunc);
        diag_file.write(reinterpret_cast<const char*>(phi_out),sizeof(double)*GS);
        diag_file.close();
        */
        for(int I=0;I<GS;I++){
            BD_out[I*Nscalars+H]=phi_x[I][0];
            BD_out[I*Nscalars+H+GS*Nscalars]=pi_x[I][0];
        }
    }

    BD_file.open("BD_file.txt", ios::out | ios::binary);
    BD_file.write(reinterpret_cast<const char*>(BD_out),sizeof(double)*GS*2*Nscalars);
    BD_file.close();
    cout<<"done!"<<endl;

    
    return 3;
}