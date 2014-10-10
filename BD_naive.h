//
//  BD_naive.h
//  BD_naive
//
//  Created by Stewart Koppell on 8/11/14.
//  Copyright (c) 2014 Stewart Koppell. All rights reserved.
//

#ifndef __BD_naive__BD_naive__
#define __BD_naive__BD_naive__

#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include "/opt/local/include/fftw3.h"
#include <ctime>

using namespace std;

double * gen_BD(int N, int Nscalars, double x_max, double H0, double * M);

#endif /* defined(__BD_naive__BD_naive__) */
