//
//  BD_naive.cpp
//  BD_naive
//
//  Created by Stewart Koppell on 8/11/14.
//  Copyright (c) 2014 Stewart Koppell. All rights reserved.
//

#include "BD_naive.h"

using namespace std;

void grand(double GRAND[], int &SEED, int rand_on)//gausian random noise
{
    SEED++;
    double PI=3.1415926535;
    srand(int(time(NULL))*SEED+SEED);
    int RANDSTRING=rand();
    double URAND=(RANDSTRING%10000)/double(10000);
    double VRAND=(((RANDSTRING-(RANDSTRING%10000+1))/10000)%10000+1)/double(10000);
    
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
    if (URAND==0 || VRAND==0) cout<<URAND<<" "<<VRAND<<" "<<endl;
    return;
}


double * gen_BD( int N, int Nscalars, double x_max, double H0, double * M){
    double PI=3.1415926535;
    double m_plank=(1e7)/3;
    int NNN=N*N*N;
    double dk=2*PI/x_max;
    int SEED=0;
    
    time_t start_time;
    time_t end_time;
    time(&start_time);

    double * BD_out;
    double * Meff;
    fftw_complex * phi_k;
    fftw_complex * phi_x;
    fftw_complex * pi_k;
    fftw_complex * pi_x;
    int * symtrack;
    int rand_on=1;
    double * phi_out;

    BD_out=new double[NNN*2*Nscalars];
    Meff=new double [Nscalars];
    for(int H=0;H<Nscalars;H++) Meff[H]=9/4.*H0+M[H];
    phi_x=new fftw_complex[NNN];
    phi_k=new fftw_complex[NNN];
    pi_x=new fftw_complex[NNN];
    pi_k=new fftw_complex[NNN];
    BD_out=new double[NNN*Nscalars*2];
    phi_out=new double[NNN];
    symtrack=new int[NNN];
    fftw_plan plan_phi;
    fftw_plan plan_pi;
    plan_phi=fftw_plan_dft_3d(N, N, N, phi_k, phi_x,FFTW_BACKWARD,FFTW_MEASURE);
    plan_pi=fftw_plan_dft_3d(N,N,N,pi_k,pi_x,FFTW_BACKWARD,FFTW_MEASURE);
    
    double kk=0;
    double rescale=1.0/sqrt(2*PI*dk*dk*dk)/sqrt(2)/m_plank;
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
        for(int I=0;I<NNN;I++) symtrack[I]=1;
        double Grand[2];
        for(int Kx=0;Kx<N;Kx++){for(int Ky=0;Ky<N;Ky++){for(int Kz=0;Kz<N;Kz++){
            kx=(N-Kx)%N; ky=(N-Ky)%N; kz=(N-Kz)%N;
            if(kx==Kx && ky==Ky && kz==Kz && symtrack[Kx+N*(Ky+N*Kz)]==1){
                symtrack[kx+N*(ky+N*kz)]=0;
                grand(Grand, SEED, rand_on);
                phi_k[Kx+N*(Ky+N*Kz)][0]=phi_k[Kx+N*(Ky+N*Kz)][0]*sqrt(Grand[0]*Grand[0]+Grand[1]*Grand[1]);
                phi_k[Kx+N*(Ky+N*Kz)][1]=0;
                grand(Grand, SEED, rand_on);
                pi_k[Kx+N*(Ky+N*Kz)][0]=phi_k[Kx+N*(Ky+N*Kz)][0]*sqrt(Grand[0]*Grand[0]+Grand[1]*Grand[1]);
                pi_k[Kx+N*(Ky+N*Kz)][1]=0;
            }
            else if(symtrack[Kx+N*(Ky+N*Kz)]==1){
                symtrack[kx+N*(ky+N*kz)]=0;
                grand(Grand, SEED, rand_on);
                phi_k[kx+N*(ky+N*kz)][1]=phi_k[kx+N*(ky+N*kz)][0]*Grand[1];
                phi_k[kx+N*(ky+N*kz)][0]=phi_k[kx+N*(ky+N*kz)][0]*Grand[0];
                grand(Grand, SEED, rand_on);
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
        for(int I=0;I<NNN;I++){
            BD_out[I+2*H*NNN]=phi_x[I][0];
            BD_out[I+2*H*NNN+NNN]=pi_x[I][0];
        }
    }

    //fstream BD_file;
    //BD_file.open("BD_file.txt", ios::out | ios::binary);
    //BD_file.write(reinterpret_cast<const char*>(BD_out),sizeof(double)*GS*2*Nscalars);
    //BD_file.close();
    

    
    time(&end_time);
    
    cout<<"BD generation done. It took "<<difftime(end_time,start_time)<<" seconds."<<endl;


    delete [] Meff;
    delete [] phi_k;
    delete [] phi_x;
    delete [] pi_k;
    delete [] pi_x;
    delete [] symtrack;
    delete [] phi_out;
    
    return BD_out;
}