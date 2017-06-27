//
// Created by chris on 5/27/2017.
//
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "IntegerPolynomialArithmetic.h"
#include "FFT.h"
using namespace std;
/**
 * Implements Rearrange on a vector of ints. Iterates 0 -> n.
 * @param x : the vector of ints.
 */
void rearrange(vector<long double> *x,int *sign,long double base){
	/*
	 * Inputs: vector of long doubles, x; valid long double base
	 * Outputs: sign of x; x itself;
	 *
	 * Assumptions: valid, non-zero base
	 */
	long double j=0; int s = 1;
	long double q; long double r;
	vector<long double>::iterator xi = x->begin();
	while(xi != x->end()){
		// j is the carry
		j = *xi + j;
		if (j < 0) s = -1; else s = 1;
		r = remainder(abs(j),base); if (r < 0) r = base + r;
		q = trunc(j / base);
		j = q;

		*xi = r;xi++;
	}
	while(j != 0){
		if (j <0) s = -1; else s = 1;
		r = remainder(abs(j),base); if (r < 0) r = base + r;
		q = trunc(j / base);
		j = q;
		x->push_back(r);
	}
	// remove trailing zeros
	if (x->size() == 0) *sign = 1;
	vector<long double>::iterator xe = x->end();
	do{
		xe--;
		if (*xe != 0) break;
		x->pop_back();
	} while(xe != x->begin());
	if (x->size() == 0) *sign = 1;
	*sign *= s;
}
/**
 * Rebases the integer polynomial to a new base.
 * @param x : integer poly to be rebased.
 * @param signX : sign of x
 * @param base : the new base of x.
 */
void rebase(vector<long double> *x,int *signX,long double old_base,long double new_base){
	if (new_base == old_base) return;
	if (new_base < old_base){
		rearrange(x,signX,new_base);
		return;
	}

	// Otherwise,
	int j = 0;long double p = 0;
	for(int i =0; i < x->size();i++){
		if ((*x)[j] > new_base) {
			j++;
		}
		p = i - j;
		(*x)[j] += powl(old_base,p)*((*x)[i]);
	}
	rearrange(x,signX,new_base);
}
/**
 * Runs addition algorithm, performs rearrangement
 * @param X : integer polynomial to be added to y
 * @param signX : sign of X
 * @param Y : integer polynomila to be added to x
 * @param signY : sign of Y
 * @param signZ : the placeholder for the sign of Z
 * @param base : the desired base of Z
 * @return Z : result of x + y
 */
vector<long double> add(vector<long double> *X,int signX,vector<long double> *Y,int signY,int *signZ,long double base){
	/*
	 * Inputs: int vector x, int vector y
	 * Outputs: x + y
	 *
	 *        Assumptions: x and y should be added; bases are correct.
	 *
	 *        Verify usages with unit tests.
	 */
	vector<long double> z (max(X->size(),Y->size()),0);
	ulong i = 0;
	long loop = min(X->size(),Y->size());
	long difference = (long) abs((double)X->size()-(double)Y->size());
	bool xgreater = X->size() > Y->size();
#pragma omp parallel for
	for(long i = 0; i < loop; i++){
		z[i] = signX*(*X)[i] + signY*(*Y)[i];
	}
	if (xgreater) {
#pragma omp parallel for
		for (long i = loop; i < difference + loop; i++) {
			z[i] = signX*(*X)[i];
		}
	} else {
#pragma omp parallel for
		for (long i = loop; i < difference + loop; i++) {
			z[i] = signY*(*Y)[i];
		}
	}
	rearrange(&z,signZ,base);
	return z;
}
/**
 * Runs subtraction algorithm via addition, performs rearrangement
 * @param X : integer polynomial to be added to y
 * @param signX : sign of X
 * @param Y : integer polynomila to be added to x
 * @param signY : sign of Y
 * @param signZ : the placeholder for the sign of Z
 * @param base : the desired base of Z
 * @return Z : result of x + y
 */
vector<long double> subtract(vector<long double> *X,int signX,vector<long double> *Y,int signY,int *signZ,long double base){
	/*
	 * Inputs: int vector x, int vector y
	 * Outputs: x - y, sign of Z
	 *
	 *
	 *        Assumptions: x and y should be added; bases are correct.
	 *
	 *        Verify usages with unit tests.
	 */
	return add(X,signX,Y,-1*signY,signZ,base);
}
/**
 * Multiplies an integer polynomial by a constant, q <= base of the polynomial.
 * Limited to int, because long double may overflow otherwise.
 * @param x : integer polynomial
 * @param signX : the sign of x
 * @param baseX : the base of x
 * @param q : the constant, q
 * @param signZ : the placeholder for the sign of z
 * @return Z: the result of qX
 */
vector<long double> multiply(vector<long double> *x,int signX,long double baseX,int q,int *signZ){
	/*
	 * Assumptions: you did not use a q value greater than the base of the polynomial.
	 *              Also, base is > 0.
	 */
	if (q == 0)	{*signZ = 1;return {};}

	// Otherwise
	vector<long double> z(x->size(),0);
	int absq = abs(q);
#pragma omp parallel for
	for(int i = 0; i < x->size(); i++){
		z[i] = (*x)[i]*absq;
	}
	int q_sign;
	if (q < 0) q_sign = -1; else q_sign = 1;
	*signZ = signX * q_sign;
	rearrange(&z,signZ,baseX);
	return z;
}
/**
 * Divides x with baseX by a number q <= baseX
 * @param x : integer polynomial
 * @param signX : sign of x
 * @param baseX : base of x
 * @param q : number to divide by; q <= base
 * @param signZ : placeholder for sign of z
 * @param r : placeholder for the remainder
 * @return quotient result; also fills in the remainder.
 */
vector<long double> divide(vector<long double> *x,int signX,long double baseX,int q,int *signZ,long double *r){
	/*
	 * Does throw a div by zero error.
	 *
	 * Assumptions: Your constant is less than or equal to your base.
	 *              Otherwise, you will need to allow for a vectored remainder...
	 *              And you should be using the polynomial division for that.
	 *
	 */
	try {if(q == 0) throw 0;}
	catch (int i){
		cout << "Error BigInt.cpp line 369: division by zero.";
	}



	if (x->size() == 0){*signZ = 1;return {};}

	vector<long double> z(x->size(),0);

	vector<long double>::const_iterator xe = x->end();
	vector<long double>::iterator ze = z.end();

	int abs_q = abs(q);
	int sign_q; if (q >= 0) sign_q = 1; else sign_q = -1;
	long double carry = 0;
	long double result;
	do{
		xe--;ze--;
		result = *xe + carry;

		(*ze) = trunc(result / abs_q);

		carry = trunc(remainder(result,abs_q));
		if (carry < 0) carry = abs_q + carry;
		carry *= baseX;

	} while(xe != x->begin());
	if (carry > 0){
		// ToDo: Does the remainder keep the sign of the final number?
		*r = carry;
	}
	// Next, the sign..
	*signZ = signX*sign_q;
	return z;
}

vector<vector<long double>> divide(vector<long double> *x,int *signX,long double *baseX,
                                   vector<long double> *y,int *signY,long double *baseY,
                                   int *signZ){
	if (*baseX != *baseY) {
		rebase(x,signX,*baseX,fmaxl(*baseX,*baseY));
		rebase(y,signY,*baseY,fmaxl(*baseX,*baseY));
		*baseX = *baseY;
	}
	if (zero(x)) {
		*signZ = 1;
		return {{},{}};
	}
	int sign_xo2;
	long double odd;
	vector<long double> x_over_2 = divide(x,*signX,*baseX,2,&sign_xo2,&odd);

	int subZ_sign;
	vector<vector<long double>> q_r = divide(&x_over_2,signX,baseX,y,signY,baseY,&subZ_sign);

	int signQ,signR;
	vector<long double> q = multiply(&(q_r[0]),sign_xo2,*baseX,2,&signQ);
	vector<long double> r = multiply(&(q_r[1]),sign_xo2,*baseX,2,&signR);

	if (odd > 1e-10 || odd < -1e-10) add(&r,signR,1,*baseX);
	if (ge_y(&r,&y)){
		vector<long double> rnew = subtract(&r,signR,y,*signY,&signR,*baseX);
		vector<long double> qnew = add(&q,signQ,1,*baseX);
		return {qnew,rnew};
	}
	return {q,r};
}