//
// Created by Chris Donlan on 5/22/2017.
//

#ifndef BIGINT3_FFT_H
#define BIGINT3_FFT_H
#include "Complex.h"
#include <elf.h>
#include<vector>

using namespace std;
vector<long double> nth_root_of_unity(long double n);
vector<long double> inverse_nth_root_of_unity(long double n);
vector<vector<long double>> FFTPrep(vector<vector<long double>> *A);
vector<vector<long double>> FFTPrep(vector<long double> *A);
void PadTo(vector<vector<long double>> *A,long double target_size);
vector<vector<long double>> PadCopy(vector<vector<long double>> *A,long double target_size);
void PadTo(vector<long double> *A,long double target_size);
vector<long double> PadCopy(vector<long double> *A,long double target_size);
void FFTMultiPrep(vector<long double> *A,vector<long double>*B);
void FFTMultiPrep(vector<vector<long double>> *A,vector<vector<long double>> *B);
vector<vector<long double>> FFT(vector<vector<long double>> *A,vector<long double> *w);
vector<vector<long double>> FFTMultiplyComplex(vector<vector<long double>> A,vector<vector<long double>> B);
vector<vector<long double>> FFTMultiplyComplex(vector<vector<long double>> *A,vector<vector<long double>> *B);
vector<long double> FFTMultiply(vector<long double> *A,vector<long double> *B);
#endif //BIGINT3_FFT_H
