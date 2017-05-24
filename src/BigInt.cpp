//
// Created by Chris Donlan on 5/16/2017.
//

#include <complex.h>
#include <math.h>
#include "BigInt.h"
#include "FFT.h"
#include "Complex.h"


using namespace std;

// ToDo: Write test: my_copy(vector<long double>::const_iterator xi, vector<long double>::const_iterator xe, vector<long double>::iterator outi);
void my_copy(vector<long double>::const_iterator xi,vector<long double>::const_iterator xe, vector<long double>::iterator outi){
	while(xi != xe){
		*outi = *xi;
		xi++;outi++;
	}
}
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

// ToDo: Write test: ddiv(long double *x, long double *y, long double *q, long double *r);
/**
 * Calculates quotient and remainder for x/y
 * @param x : numerator pointer
 * @param y : divisor pointer
 * @param q : quotient pointer
 * @param r : remainder pointer
 */
inline void ddiv(long double *x,long double *y,long double *q,long double *r){
	*q = *x / *y;
	*r = modfl(*x,y);
	int xx = 0;
}
// ToDo: Write test: ddiv(long double *x,int *y, long double *q, long double *r);
/**
 * Calculates quotient and remainder for x/y
 * @param x : numerator pointer
 * @param y : divisor pointer
 * @param q : quotient pointer
 * @param r : remainder pointer
 */
inline void ddiv(long double *x,int *y,long double *q,long double *r){
	*q = *x / *y;
	*r = trunc(abs(remainder(*x,*y)));
}
// ToDo: Write test: ddiv(vector<long double>::iterator *xi, long double *y, long double *q, long double *r);
/**
 * Calculates quotient and remainder for *xi / y
 * @param xi : iterator whose value pointed to is divided by y
 * @param y  : divisor pointer
 * @param q  : quotient pointer
 * @param r  : remainder pointer
 */
inline void ddiv(vector<long double>::iterator * xi,long double *y,long double *q,long double *r){
	//todo: optimize this for one function call!
	*q = trunc(**xi / *y);
	*r = trunc(abs(remainder(**xi,*y)));
}

// ToDo: Write test: ddiv(vector<long double>::const_iterator * xi, long double *y, long double *q,long double *r);
/**
 * Calculates quotient and remainder for *xi / y
 * @param xi : const_iterator whose value pointed to is divided by y
 * @param y  : divisor pointer
 * @param q  : quotient pointer
 * @param r  : remainder pointer
 */
inline void ddiv(vector<long double>::const_iterator * xi,long double *y,long double *q,long double *r){
	//todo: optimize this for one function call!
	*q = trunc(**xi / *y);
	*r = trunc(abs(remainder(**xi,*y)));
}
// ToDo: Write test: ddiv(vector<long double>::const_iterator *xi, int *y, long double *q, long double *r);
/**
 * Calculates quotient and remainder for *xi / y
 * @param xi : const_iterator whose value pointed to is divided by y
 * @param y  : divisor pointer
 * @param q  : quotient pointer
 * @param r  : remainder pointer
 */
inline void ddiv(vector<long double>::const_iterator * xi,int *y,long double *q,long double *r){
	//todo: optimize this for one function call!
	*q = trunc(**xi / *y);
	*r = trunc(abs(remainder(**xi,*y)));
}
// ToDo: Write Test: rearrange_subproblem(vector<long double>::iterator *xi,long double *B,long double * carry)
/**
 * Insertion point for the rearrange function.
 * @param xi : the iterator used in the rearrange
 * @param B  : the base of the BigInt
 * @param carry : the pointer to the value to be carried through the rearrange
 * @return : sign of the subproblem (+|-);
 */
int  rearrange_subproblem(vector<long double>::iterator *xi,long double *B,long double * carry){
	/*
	// todo: make sure I choose a base such that I will not get overflow.
	// todo: possibly, optimize
	//
	// Input: BigInt vector<long double> iterator xi,
	//        BigInt Base pointer *B,
	//        The "carry" placeholder pointer
	//
	// Output: updates carry, and outputs the sign of number at position xi after adding carry.
	// Todo: handle the zero case--can't have a negative if the number is zero.
	// ----------------------------------------------------------------------------------------*/
	int sign;
	if (**xi < 0) sign = -1; else sign = 1;
	**xi += (long double)((long long)*carry);


	long double xiv = abs(**xi);  //
	if (abs(**xi) > *B) {
		long double r;
		ddiv(xi,B,carry,&r);
		*carry *= sign;
		**xi = (long double)(trunc(abs(r)));
	} else if (abs(**xi) == *B){
		**xi = 0;
		*carry = sign;

	} else *carry = 0;
	return sign;
}

// ToDo: Write Test: rearrange(BigInt *X);
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

		// ToDo: turn this into one action.
		r = remainder(j,X->base);
		q = trunc(j / X->base);
		if (q == 0) r = j;
		j = q;

		*xi = r;xi++;
	}
	while(j != 0){
		if (j <0) sign = -1; else sign = 1;
		r = remainder(j,X->base);
		q = trunc(j / X->base);
		if (q == 0) r = j;
		j = q;

		X->number.push_back(r);
	}
	X->sign *= sign;
	clean_trailing_zeros(X);
	if (X->number.size() == 0) X->sign = 1;

}
// ToDo: Write Test: rearrange(BigInt *X,int (*subproblem));
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
BigInt SimpleArithmetic(BigInt *X,BigInt *Y, long double (*subproblem)(int *signX,vector<long double>::const_iterator *xiter, int *signY,vector<long double>::const_iterator *yiter)){
	BigInt Z(max(X->number.size(),Y->number.size()),X->base,1);
	ulong loop_size = min(X->number.size(),Y->number.size());
	ulong i = 0;
	vector<long double>::const_iterator xi = X->number.begin();
	vector<long double>::const_iterator yi = Y->number.begin();
	vector<long double>::iterator zi = Z.number.begin();

	// ToDo: Add: may want to make a parallel version
	while(i < loop_size){
		*zi = subproblem(&(X->sign),&xi,&(Y->sign),&yi);

		zi++;xi++;yi++;i++;
	}
	// ToDo: Add: make sure the copying works properly.
	my_copy(xi,X->number.end(),Z.number.begin()+loop_size);
	my_copy(yi,Y->number.end(),Z.number.begin()+loop_size);

	// Rearrange update's Z's sign.
	rearrange(&Z);
	// delete trailing zeros
	// ...if all zeros, set sign to 1
	clean_trailing_zeros(&Z);
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
inline long double AddSubproblem(int *signX,vector<long double>::const_iterator *xiter, int *signY,vector<long double>::const_iterator *yiter){
	return *signX * (**xiter) + *signY*(**yiter);
}
/**
 * Add subproblem for use in the function "SimpleArithmetic"
 * @param signX : the "sign" member of the BigInt X (left)
 * @param xiter : the iterator for BigInt X (left BigInt)
 * @param signY : the "sign" member of the BigInt Y (right)
 * @param yiter : the iterator for BigInt Y (right BigInt)
 * @return : BigInt Z = X - Y read left to right.
 */
inline long double SubtractSubproblem(int *signX,vector<long double>::const_iterator *xiter, int *signY,vector<long double>::const_iterator *yiter){
	return (*signX * (**xiter)) - (*signY*(**yiter));
}

void rebase(BigInt *A, BigInt *B){
	if (A->base != B->base){

		//ToDo: Rebase!!!
		long double old_base = B->base;
		long double a_base = A->base;

		if (a_base < old_base) {
			B->base = A->base;
			B->rearranged = false;
			B->update();
		}
		else {
			// ToDo implement this.
		}
	}
}
// Todo:Test: Add(BigInt *X,BigInt *Y)
// Todo:Runtime and Complexity Test: Add(BigInt *X,BigInt *Y)
/**
 * Evaluates X + Y = Z, left to right, in an immutable fashion (creates new BigInt, Z)
 * @param X : Left hand BigInt Pointer
 * @param Y : Right hand BigInt Pointer
 * @return Z : solution to X + Y
 */
BigInt Add(BigInt *X,BigInt *Y){
	rebase(X, Y);
	return SimpleArithmetic(X,Y,AddSubproblem);
}

// Todo:Test: Subtract(BigInt *X, BigInt *Y)
// ToDo:Runtime and Complexity Test: Subtract(BigInt *X, BigInt *Y)
/**
 * Evaluates X - Y = Z, left to right, in an immutable fashion (creates new BigInt, Z)
 * @param X : Left hand BigInt Pointer
 * @param Y : Right hand BigInt Pointer
 * @return Z : solution to X - Y, read left to right.
 */
BigInt Subtract(BigInt *X,BigInt *Y){
	rebase(X, Y);
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

// ToDo:Document: MultiplySubproblem(...)
inline long double MultiplySubproblem(int *signX,vector<long double>::const_iterator *xiter,
                                      int *signY,vector<long double>::const_iterator *yiter){
	return (*signX * (**xiter)) * (*signY*(**yiter));
}

// ToDo: Document: Multiply(BigInt *X,int q)
// Todo:Test: Multiply(BigInt *X, int q)
// Todo: Runtime and Complexity Test: Multiply(BigInt *X, int q);
BigInt Multiply(BigInt *X,int q){
	if (q == 0){
		vector<long double> n;
		BigInt Z(n,X->base,1);
		return Z;
	}
	vector<long double> n(X->number.size(),0);
	BigInt Z(n,X->base,1);


	int i;
#pragma omp parallel for
	for(i = 0; i < X->number.size(); i++){
		Z.number[i] = X->number[i]*q;
	}

	int q_sign;
	if (q < 0) q_sign = -1; else q_sign = 1;
	Z.sign = X->sign * q_sign;
	Z.rearranged = false;
	Z.update();
	return Z;
}

// ToDo: Document: Divide(BigInt *X, int q, long double * remainder)
// Todo:Test: Divide(BigInt*X, int q, long double * remainder)
// ToDo:Runtime and Complexity Test: Divide(BigInt*X, int q, long double * remainder)
BigInt Divide(BigInt *X, int q,long double * remainder){
	BigInt Z(X->number.size(),X->base,1);
	vector<long double>::const_iterator xe = X->number.end();
	vector<long double>::iterator ze = Z.number.end();

	int abs_q = abs(q);
	int sign_q; if (q >= 0) sign_q = 1; else sign_q = -1;
	long double carry = 0;
	long double result;
	do{
		xe--;ze--;
		result = *xe + carry;
		ddiv(&result,&q,&(*ze),&carry);
		carry *= X->base;

	} while(xe != X->number.begin());
	if (carry > 0){
		// ToDo: Does the remainder keep the sign of the final number?
		*remainder = carry;
	}
	// Next, the sign..
	Z.sign *= sign_q;
}

// ToDo: Document: Multiply
// ToDo:Test: Multiply(BigInt *A, BigInt *B)
// ToDo:Runtime and Complexity test: Multiply(BigInt *A, BigInt *B)
BigInt Multiply(BigInt *A,BigInt *B){
	// ToDo: finish this.
	// ToDo: rebase to optimal base.


	
	// ToDo: both A and B need to have the same base.
	// for now assume that A and B have a good base, and that the bases are equivalent.
	vector<vector<long double>> Ac = Real2Complex(&(A->number));
	vector<vector<long double>> Bc = Real2Complex(&(B->number));
	vector<vector<long double>> Cc = FFTMultiplyComplex(Ac,Bc);

	vector<long double> C  = Complex2RealPolar(&Cc);

	BigInt Out(C,A->base,A->sign*B->sign);
	return Out;
}





// ToDo: Barret's algorithm
//  ... cant go further before I lock down what I have done so far.



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
