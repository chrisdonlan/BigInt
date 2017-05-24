//
// Created by Chris Donlan on 5/22/2017.
//

#ifndef BIGINT3_COMPLEX_H
#define BIGINT3_COMPLEX_H
#include <vector>
#define M_PIl  3.141592653589793238462643383279502884L
using namespace std;
vector<long double> polar(vector<long double> *c); // Done
vector<vector<long double>> polar(vector<vector<long double>> *c); // Done
vector<long double> cartesian(vector<long double> *polar_complex); // Done
vector<vector<long double>> cartesian(vector<vector<long double>> *p); // Done
vector<long double> subtract_cartesian(vector<long double>*a,vector<long double>*b); // Done
vector<long double> add_cartesian(vector<long double>*a,vector<long double>*b); // Done
vector<long double> multiply_polar(vector<long double> *a,vector<long double> *b); // Done

vector<long double> multiply_cartesian(vector<long double> *a,vector<long double> *b); // Done
vector<long double> subtract_polar(vector<long double> *a,vector<long double> *b); // Done
vector<long double> add_polar(vector<long double>*a,vector<long double> *b); // Done
vector<vector<long double>> multiply_polar(vector<vector<long double>> *A, vector<vector<long double>> *B); // Done
vector<vector<long double>> multiply_cartesian(vector<vector<long double>> *A, vector<vector<long double>> *B); // Done
vector<long double> multiply_polar(long double *a,vector<long double> *b);  // Done
vector<long double> multiply_cartesian(long double *a,vector<long double> *b);// Done

vector<vector<long double>> multiply_polar(long double *a,vector<vector<long double>> *b);// Done
vector<vector<long double>> multiply_cartesian(long double *a,vector<vector<long double>> *b); // Done

vector<vector<long double>> multiply_polar(vector<long double> *a,vector<vector<long double>> *b);// Done
vector<vector<long double>> multiply_cartesian(vector<long double> *a,vector<vector<long double>> *b); // Done

vector<vector<long double>> Real2Complex(vector<long double> *A); // Done
vector<long double> Complex2RealCartesian(vector<vector<long double>> *A); // Done
vector<long double> Complex2RealPolar(vector<vector<long double>> *A); // Done
#endif //BIGINT3_COMPLEX_H
