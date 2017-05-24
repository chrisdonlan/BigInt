//
// Created by Chris Donlan on 5/22/2017.
//

#ifndef BIGINT3_FFT_H
#define BIGINT3_FFT_H
#include "Complex.h"
#include <elf.h>
#include<vector>

using namespace std;
vector<vector<long double>> FFT(vector<vector<long double>> *A,vector<long double> *w);
vector<vector<long double>> FFTMultiplyComplex(vector<vector<long double>> A,vector<vector<long double>> B);
#endif //BIGINT3_FFT_H
