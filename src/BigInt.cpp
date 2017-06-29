//
// Created by Chris Donlan on 5/16/2017. copyright, all rights reserved
//

#include <complex.h>
#include <math.h>
#include <cmath>
#include "BigInt.h"
#include "FFT.h"
#include "Complex.h"


using namespace std;

void clean_trailing_zeros(BigInt* Z){
	if (Z->number.size() > 0) {
		vector<long double>::const_iterator ze = Z->number.end();
		bool non_empty_element_detected = false;
		while (!non_empty_element_detected) {
			ze--;
			if (*ze != 0) non_empty_element_detected = true;
			if (!non_empty_element_detected) Z->number.pop_back();
		}
	}
}
/**
 * Implements Rearrange on a BigInt. Iterates 0 -> n.
 * @param X : the BigInt Pointer
 */
void rearrange(BigInt *X){
	try {
		if (X->base == 0) throw 0;
		if (X->base == 1) throw 1;
	}
	catch (int e){
		if (e == 0)
			cout << "Singularity in BigInt: cannot rebase to 0.";
		if (e == 1)
			cout << "Singularity in BigInt: cannot rebase to 1.";
	}
	long double j=0; int sign = 1;
	long double q; long double r;
	vector<long double>::iterator xi = X->number.begin();
	while(xi != X->number.end()){
		// j is the carry
		j = *xi + j;
		if (j < 0) sign = -1; else sign = 1;
		r = remainder(abs(j),X->base); if (r < 0) r = X->base + r;
		q = trunc(j / X->base);
		j = q;

		*xi = r;xi++;
	}
	while(j != 0){
		if (j <0) sign = -1; else sign = 1;
		r = remainder(abs(j),X->base); if (r < 0) r = X->base + r;
		q = trunc(j / X->base);
		j = q;

		X->number.push_back(r);
	}
	X->sign *= sign;
	clean_trailing_zeros(X);
	if (X->number.size() == 0) X->sign = 1;
}
/**
 * Implements rearrange on a BigInt. Iterates 0 -> n.
 * @param X : BigInt pointer
 * @param subproblem : function to implement at each step along the BigInt
 */
void rearrange(BigInt *X,int (*subproblem)(vector<long double>::iterator *xi,long double *B,long double * carry)){
	long double carry = 0; int sign = 1; bool all_zeros = true;
	vector<long double>::iterator xi = X->number.begin();
	while(xi != X->number.end()){
		// carry carries through the sign. What matters is the sign of the last carry.
		sign *= subproblem(&xi,&(X->base),&carry);
		if (all_zeros && *xi > 0) all_zeros = false;
		xi++;
	}
	if (carry != 0){
		// No need to update "all_zeros" since I return from the function here.
		if (carry < 0) sign *= -1;
		X->number.push_back(trunc((int)abs(carry)));

		// ToDo: make sure all of these signs are correct.
		X->sign *= sign;
		return;
	}
	if (all_zeros) X->sign = 1;
}
//region Add and Subtract (Simple Arithmetic)
/**
 * Baseline function for BigInt Add and Subtract
 * @param X : BigInt pointer to the left number, reading (X +|- Y) left to right
 * @param Y : BigInt pointer to the right number, reading (X +|- Y) left to right
 * @param subproblem : The actual addition or subtraction step (AddSubproblem|SubtractSubproblem)
 * @return : the BigInt that is the value of the operation
 */
BigInt SimpleArithmetic(BigInt *X,BigInt *Y, long double (*subproblem)(long double x, long double y)){
	vector<long double> z (max(X->number.size(),Y->number.size()),0);
	ulong loop_size = min(X->number.size(),Y->number.size());
	ulong i = 0;

	// loop size
	long loop = min(X->number.size(),Y->number.size());
	// difference
	long difference = abs((long) (X->number.size() - Y->number.size()));
	bool xgreater = X->number.size() > Y->number.size();
#pragma omp parallel for
	for(long i = 0; i < loop; i++){
		z[i] = subproblem(X->sign*(X->number[i]),Y->sign*(Y->number[i]));
	}
	if (xgreater) {
#pragma omp parallel for
		for (long i = loop; i < difference + loop; i++) {
			z[i] = X->number[i] * (X->sign);
		}
	} else {
#pragma omp parallel for
		for (long i = loop; i < difference + loop; i++) {
			z[i] = Y->number[i] * (Y->sign);
		}
	}
	BigInt Z(z,X->base,1);
	return Z;
}
/**
 * Add subproblem for use in the function "SimpleArithmetic"
 * @param signX : the "sign" member of the BigInt X (left)
 * @param xiter : the iterator for BigInt X (left BigInt)
 * @param signY : the "sign" member of the BigInt Y (right)
 * @param yiter : the iterator for BigInt Y (right BigInt)
 * @return : BigInt Z = X + Y read left to right.
 */
inline long double AddSubproblem(long double x, long double y){
	return x + y;
}
/**
 * Add subproblem for use in the function "SimpleArithmetic"
 * @param signX : the "sign" member of the BigInt X (left)
 * @param xiter : the iterator for BigInt X (left BigInt)
 * @param signY : the "sign" member of the BigInt Y (right)
 * @param yiter : the iterator for BigInt Y (right BigInt)
 * @return : BigInt Z = X - Y read left to right.
 */
inline long double SubtractSubproblem(long double x, long double y){
	return x - y;
}

void rebase(BigInt *A, long double base){
	if (base <= A->base){
		A->base = base;
		A->update();
		return;
	} else {
		long double old_base = A->base;
		vector<long double> new_number(1,0);
		int j = 0;long double p = 0;
		for(int i =0; i < A->number.size();i++){
			if (new_number[j] < base) {
				j++;
				new_number.push_back(0);
			}
			p = i - j;
			new_number[j] += powl(old_base,p)*(A->number[i]);
		}
		A->number = new_number;
		A->update();
		return;
	}
	// so, I start accumulating numbers in A until it is >= base, or numbers run out!
}
// Todo:Runtime and Complexity Test: Add(BigInt *X,BigInt *Y)
/**
 * Evaluates X + Y = Z, left to right, in an immutable fashion (creates new BigInt, Z)
 * @param X : Left hand BigInt Pointer
 * @param Y : Right hand BigInt Pointer
 * @return Z : solution to X + Y
 */
BigInt Add(BigInt *X,BigInt *Y){
	// ToDo: parallelize function calls!
	rebase(X,max(X->base,Y->base));
	rebase(Y,max(X->base,Y->base));
	return SimpleArithmetic(X,Y,AddSubproblem);
}
// ToDo:Runtime and Complexity Test: Subtract(BigInt *X, BigInt *Y)
/**
 * Evaluates X - Y = Z, left to right, in an immutable fashion (creates new BigInt, Z)
 * @param X : Left hand BigInt Pointer
 * @param Y : Right hand BigInt Pointer
 * @return Z : solution to X - Y, read left to right.
 */
BigInt Subtract(BigInt *X,BigInt *Y){
	rebase(X,max(X->base,Y->base));
	rebase(Y,max(X->base,Y->base));
	return SimpleArithmetic(X,Y,SubtractSubproblem);
}
//endregion
/* Original Add...
BigInt Add(BigInt *X,BigInt *Y){
	BigInt Z(max(X->number.size(),Y->number.size()),X->base,1);
	ulong loop_size = min(X->number.size(),Y->number.size());
	ulong i = 0;
	vector<long double>::const_iterator xi = X->number.begin();
	vector<long double>::const_iterator yi = Y->number.begin();
	vector<long double>::iterator zi = Z.number.begin();


	while(i < loop_size){
		*zi = X->sign *  *xi + Y->sign * *yi;
		zi++;xi++;yi++;i++;
	}

	my_copy(xi,X->number.end(),Z.number.begin()+loop_size);
	my_copy(yi,Y->number.end(),Z.number.begin()+loop_size);

	// Rearrange update's Z's sign.
	rearrange(&Z);
	// delete trailing zeros
	// ...if all zeros, set sign to 1
	clean_trailing_zeros(&Z);
	return Z;
}
*/
/* Original Subtract
BigInt Subtract(BigInt *X,BigInt *Y){
	BigInt Z(max(X->number.size(),Y->number.size()),X->base,1);
	ulong loop_size = min(X->number.size(),Y->number.size());
	ulong i = 0;
	vector<long double>::const_iterator xi = X->number.begin();
	vector<long double>::const_iterator yi = Y->number.begin();
	vector<long double>::iterator zi = Z.number.begin();

	// ToDo: Subtract: may want to make a parallel version
	while(i < loop_size){
		*zi = X->sign *  *xi - Y->sign * *yi;
		zi++;xi++;yi++;i++;
	}
	// ToDo: Subtract: make sure the copying works properly.
	my_copy(xi,X->number.end(),Z.number.begin()+loop_size);
	my_copy(yi,Y->number.end(),Z.number.begin()+loop_size);

	// Rearrange update's Z's sign.
	rearrange(&Z);
	// delete trailing zeros
	// ...if all zeros, set sign to 1
	clean_trailing_zeros(&Z);
	return Z;
}
*/

// ToDo: Document: Multiply(BigInt *X,int q)
// Todo: Runtime and Complexity Test: Multiply(BigInt *X, int q);
BigInt Multiply(BigInt *X,int q){
	if (q == 0){
		vector<long double> n;
		BigInt Z(n,X->base,1);
		return Z;
	}
	vector<long double> z(X->number.size(),0);
	int absq = abs(q);
#pragma omp parallel for
	for(int i = 0; i < X->number.size(); i++){
		z[i] = X->number[i]*absq;
	}
	BigInt Z(z,X->base,1);
	int q_sign;
	if (q < 0) q_sign = -1; else q_sign = 1;
	Z.sign = X->sign * q_sign;
	return Z;
}

// ToDo: Document: Divide(BigInt *X, int q, long double * remainder)
// ToDo:Runtime and Complexity Test: Divide(BigInt*X, int q, long double * remainder)
BigInt Divide(BigInt *X, int q,long double * r){
	try {if(q == 0) throw 0;}
	catch (int i){
		cout << "Error BigInt.cpp line 369: division by zero.";
	}
	if (X->number.size() == 0){BigInt Z({0},X->base,1); return Z;}

	vector<long double> z(X->number.size(),0);

	vector<long double>::const_iterator xe = X->number.end();
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
		carry *= X->base;

	} while(xe != X->number.begin());
	if (carry > 0){
		// ToDo: Does the remainder keep the sign of the final number?
		*r = carry;
	}
	// Next, the sign..
	BigInt Z(z,X->base,X->sign);
	Z.sign *= sign_q;
	return Z;
}

// ToDo:Runtime and Complexity test: Multiply(BigInt *A, BigInt *B)
/// Multiply uses FFT to accomplish the multiply.  Assumes that the bases of A and B have an optimal base,
/// and sufficiently small to guarantee that integer multiplication will be correct.
/// rule of thumb is that INT32_MAX / 256 is a good base..Will have to verify this, though.
/// accuracy changes with O(B^2 * precision of number type * length of number vector);
/// \param A : A BigInt
/// \param B : B BigInt
/// \return Z, the result of A*B
BigInt Multiply(BigInt *A,BigInt *B){
	// for now assume that A and B have a good base, and that the bases are equivalent.
	vector<vector<long double>> Ac = Real2Complex(&(A->number));
	vector<vector<long double>> Bc = Real2Complex(&(B->number));
	vector<long double> C = FFTMultiply(&(A->number),&(B->number));
	// round to nearest integer
#pragma omp parallel for
	for(int i = 0; i < C.size(); i++) C[i] = round(C[i]);

	BigInt Out(C,A->base,A->sign*B->sign);
	return Out;
}


bool zero(vector<long double> *x){
	if (x->size() == 0) return true;
	volatile bool zero = true;
#pragma omp parallel for shared(zero)
	for(int i = 0; i < x-> size(); i++){
		if (zero) continue;
		if ((*x)[i] != 0) zero = false;
	}
	return zero;
}
bool ge_one(vector<long double> *x){
	if (x->size() == 0) return false;
	if ((*x)[0] >= 1) return true;

	volatile bool geone = false;
#pragma omp parallel for shared(geone)
	for(int i=0; i < x->size(); i++){

	}
}
vector<vector<long double>> divide(vector<long double> *x,vector<long double> *y){
	/*
	 * Input: Two n*sizeof(long double) bit numbers x, y, y >= 1;
	 * Output: quotient, remainder of x divided by y
	 */
	if (zero(x)) return {{0},{0}};
	vector<long double> x_over_2 = divide(x,2);
	vector<vector<long double>> q_r = divide(&x_over_2,y);

	vector<long double> q = multiply(&(q_r[0]),2);
	vector<long double> r = multiply(&(q_r[1]),2);
	if (remainder(trunc((*x)[0]),2)!= 0) r = add(r,1);
	if (ge_y(&r,&y)){
		r = subtract(&r,y);
		q = add(q,1);
	}
	return {q,r};
}
// ToDo: Barret's algorithm
//  ... cant go further before I lock down what I have done so far.
BigInt Divide(BigInt *A,BigInt *B){
	vector<long double> a(A->number.size());
	copy(A->number.begin(),A->number.end(),a.begin());
	vector<long double> b(B->number.size());
	copy(B->number.begin(),B->number.end(),b.begin());
	try{
		if(not ge_one(&b))
			if(zero(&b)) throw 0;
			else throw 1;
	} catch (int i){
		if (i == 0) cout<< "divisor is zero. must be greater than or equal to one.";
		if (i == 1) cout<< "divisor is  <1.  Must be >= 1.";
	}

}


void BigInt::update(){
	if(!(this->rearranged)) {
		rearrange(this);
		this->rearranged = true;
	}
}
BigInt::BigInt(ulong size,long double B,int s){
	vector<long double> n (size);
	number = n;
	base = B;
	sign = s;
	rearranged = false;
	this->update();
}
BigInt::BigInt(vector<long double> num,long double B,int s){
	vector<long double> n(num.size());
	copy(num.begin(),num.end(),n.begin());
	number = n;
	base = B;
	sign = s;
	rearranged = false;
	// rebase the vector
	this->update();
}
BigInt::BigInt(){
	vector<long double> n;
	number = n;
	base = INT32_MAX;
	sign = 1;
	rearranged = false;
	this->update();
}
bool vector_equality(vector<long double> * x,vector<long double> * y){
	if (x->size() != y->size()) return false;

	vector<long double>::const_iterator xi = x->begin();
	vector<long double>::const_iterator yi = y->begin();
	while(xi != x->end()){ if(*xi != *yi) return false; xi++;yi++;}
	return true;
}
bool BigInt::equals(BigInt *X){
	this->update(); X->update();
	return X->sign == this->sign && vector_equality(&(this->number),&(X->number));
}
bool BigInt::equals(BigInt X){
	this->update();X.update();
	return X.sign == this->sign && vector_equality(&(this->number),&(X.number));
}
bool BigInt::equals(vector<long double> *a){
	this->update();
	BigInt X(*a,this->base,1);
	X.update();
	return vector_equality(&(this->number),&(X.number));
}
bool BigInt::equals(vector<long double> a){
	this->update();
	BigInt X(a,this->base,1);
	X.update();
	return vector_equality(&(this->number),&(X.number));
}
