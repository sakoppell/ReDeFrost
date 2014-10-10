//
//  ReDeFrost.h
//  
//
//  Created by Stewart Koppell on 5/11/14.
//
//

#ifndef ____ReDeFrost__
#define ____ReDeFrost__

#include <algorithm>
#include <complex>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include "BD_naive.h"

using namespace std;

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


#endif /* defined(____ReDeFrost__) */
