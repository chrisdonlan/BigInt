//
// Created by Chris Donlan on 5/22/2017.
//

#include <iostream>
#include <math.h>
#include "Complex.h"
vector<long double> polar(vector<long double> *c){
	long double r = sqrtl(powl((*c)[0],2)+powl((*c)[1],2));
	long double theta = atanl((*c)[1]/(*c)[0]);
	// note that the range of arctan is -pi/2 -> pi/2

	bool x = (*c)[0] >=0,i = (*c)[1] >=0;
	// x = true, c= true, first quadrant
	if      (x&&i) return {r,          theta}; // Q1
	else if (x)    return {r,2*M_PIl + theta}; // Q4
	else if (i)    return {r,  M_PIl + theta}; // Q2
	else           return {r,  M_PIl + theta}; // Q3
}
vector<vector<long double>> polar(vector<vector<long double>> *c){
	vector<vector<long double>> p(c->size());
#pragma omp parallel for
	for(int i = 0; i < c->size(); i++)
		p[i] = polar(&(*c)[i]);
	return p;
}
vector<long double> cartesian(vector<long double> *polar_complex){
	return {(*polar_complex)[0]*cosl((*polar_complex)[1]),(*polar_complex)[0]*sinl((*polar_complex)[1])};
}
vector<vector<long double>> cartesian(vector<vector<long double>> *p){
	vector<vector<long double>> c(p->size());
#pragma omp parallel for
	for(int i = 0; i < p->size(); i++){
		c[i] = cartesian(&(*p)[i]);
	}
	return c;
}
vector<long double> subtract_cartesian(vector<long double>*a,vector<long double>*b){
	return {(*a)[0] - (*b)[0],(*a)[1] - (*b)[1]};
}
vector<long double> subtract_polar(vector<long double> *a,vector<long double> *b){
	vector<long double> ca = cartesian(a);
	vector<long double> cb = cartesian(b);
	vector<long double> cc = subtract_cartesian(&ca,&cb);
	return polar(&cc);
}

vector<long double> add_cartesian(vector<long double>*a,vector<long double>*b){
	return {(*a)[0] + (*b)[0],(*a)[1] + (*b)[1]};
}
vector<long double> add_polar(vector<long double>*a,vector<long double> *b){
	vector<long double> ca = cartesian(a);
	vector<long double> cb = cartesian(b);
	vector<long double> cc = add_cartesian(&ca,&cb);
	return polar(&cc);
}


vector<long double> multiply_polar    (vector<long double> *a,vector<long double> *b){
	return {(*a)[0]*(*b)[0],(*a)[1] + (*b)[1]};
}
vector<long double> multiply_cartesian(vector<long double> *a, vector<long double> *b){
	vector<long double> pa = polar(a),pb = polar(b);
	vector<long double> pc = multiply_polar(&pa,&pb);
	return cartesian(&pc);
}

vector<long double> multiply_polar    (long double *a, vector<long double> *b){
	return {(*b)[0]**a,(*b)[1]};
}
vector<long double> multiply_cartesian(long double *a, vector<long double> *b){
	return {(*b)[0]**a,(*b)[1]**a};
}
// ToDo: come up with a generic function for these two repetitive functions: multiply_polar/cartesian(long double a,vector<vector<long double>> *b);
vector<vector<long double>> multiply_polar    (long double *a,vector<vector<long double>> *b){
	vector<vector<long double>> c(b->size());
#pragma omp parallel for
	for(int i = 0; i < b->size(); i++)
		c[i] = multiply_polar(a,&(*b)[i]);
	return c;
}
vector<vector<long double>> multiply_cartesian(long double *a,vector<vector<long double>> *b){
	vector<vector<long double>> c(b->size());
#pragma omp parallel for
	for(int i = 0; i < b-> size(); i++)
		c[i] = multiply_cartesian(a,&(*b)[i]);
	return c;
}

vector<vector<long double>> multiply_polar    (vector<long double> *a, vector<vector<long double>> *b) {
	try{
		if (a->size() != b->size()) throw 0;
	} catch (int i){ cout << "Error: scalar vector must be the same size as the complex vector to multiply.";}
	vector<vector<long double>> c(b->size());
#pragma omp parallel for
	for(int i = 0; i < a->size(); i++)
		c[i] = multiply_polar(&(*a)[i],&(*b)[i]);
	return c;
}
vector<vector<long double>> multiply_cartesian(vector<long double> *a,vector<vector<long double>>*b) {
	try{
		if (a->size() != b->size()) throw 0;
	} catch (int i) { cout << "Error: scalar vector must be the same size as the complex vector to multiply.";}

	vector<vector<long double>> c(b->size());
#pragma omp parallel for
	for(int i =0; i < a-> size(); i++)
		c[i] = multiply_cartesian(&(*a)[i],&(*b)[i]);
	return c;
}



vector<vector<long double>> multiply_polar    (vector<vector<long double>> *A, vector<vector<long double>> *B) {
	try{if(A->size() != B->size()) throw 0;}
	catch (int e) { cout<<"Error, Complex: Hadamard Product: Input vectors must be same size.";}

	vector<vector<long double>> Out(A->size());
	int i;
#pragma omp parallel for
	for(i = 0; i < A->size(); i++){Out[i] = multiply_polar(&((*A)[i]),&((*B)[i]));}

	return Out;
}
vector<vector<long double>> multiply_cartesian(vector<vector<long double>> *A, vector<vector<long double>> *B) {
	try{if(A->size() != B->size()) throw 0;}
	catch (int e) { cout<<"Error, Input vectors must be same size.";}

	vector<vector<long double>> Out(A->size());
	int i;
#pragma omp parallel for
	for(i = 0; i < A->size(); i++){
		vector<long double> pA = polar(&((*A)[i]));
		vector<long double> pB = polar(&((*B)[i]));
		vector<long double> pC = multiply_polar(&pA,&pB);
		Out[i] = cartesian(&pC);
	}

	return Out;
}


vector<long double> Complex2RealCartesian(vector<vector<long double>> *A) {
	vector<long double> B(A->size());
	int i;
#pragma omp parallel for
	for (i = 0; i < A->size(); i++)B[i] = (*A)[i][0];
	return B;
}
vector<long double> Complex2RealPolar(vector<vector<long double>> *A) {
	vector<long double> B(A->size());
	int i;
#pragma omp parallel for
	for (i = 0; i < A->size(); i++) B[i] = (*A)[i][0]* cosl((*A)[i][1]);
	return B;
}
vector<vector<long double>> Real2Complex(vector<long double> *A){
	vector<vector<long double>> B(A->size());
	int i;
#pragma omp parallel for
	for(i = 0; i < A->size(); i++){B[i] = {(*A)[i],0};}
	return B;
}

