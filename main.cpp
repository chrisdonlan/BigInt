#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>

#define M_PIl  3.141592653589793238462643383279502884L

using namespace std;
void my_copy(vector<long double>::const_iterator xi,vector<long double>::const_iterator xe, vector<long double>::iterator outi){
	while(xi != xe){
		*outi = *xi;
		xi++;outi++;
	}
}

struct BigInt
{

	BigInt(ulong size,long double B,int s){
		vector<long double> n (size);
		number = n;
		base = B;
		sign = s;
	}
	BigInt(vector<long double> num,long double B,int s){
		vector<long double> n(num.size());
		copy(num.begin(),num.end(),n.begin());
		number = n;
		base = B;
		sign = s;
	}
	vector<long double> number;
	long double base;
	int sign;
};
void clean_trailing_zeros(BigInt* Z){
	vector<long double>::const_iterator ze = Z->number.end();
	bool non_empty_element_detected = false;
	while(!non_empty_element_detected){
		ze--;
		if (*ze != 0) non_empty_element_detected = true;
		if (!non_empty_element_detected) Z->number.pop_back();
	}
}
/**
 * Calculates quotient and remainder for x/y
 * @param x : numerator pointer
 * @param y : divisor pointer
 * @param q : quotient pointer
 * @param r : remainder pointer
 */
inline void ddiv(long double *x,long double *y,long double *q,long double *r){
	*q = *x / *y;
	*r = remainder(*x,*y);
}
/**
 * Calculates quotient and remainder for x/y
 * @param x : numerator pointer
 * @param y : divisor pointer
 * @param q : quotient pointer
 * @param r : remainder pointer
 */
inline void ddiv(long double *x,int *y,long double *q,long double *r){
	*q = *x / *y;
	*r = remainder(*x,*y);
}
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
	*r = remainder(**xi,*y);
}
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
	*r = remainder(**xi,*y);
}
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
	*r = remainder(**xi,*y);
}
/**
 * Insertion point for the rearrange function.
 * @param xi : the iterator used in the rearrange
 * @param B  : the base of the BigInt
 * @param carry : the pointer to the value to be carried through the rearrange
 * @return : sign of the subproblem (+|-);
 */
int  rearrange_subproblem(vector<long double>::iterator *xi,long double *B,long double * carry){
	// todo: make sure I choose a base such that I will not get overflow.
	// todo: possibly, optimize
	//
	// Input: BigInt vector<long double> iterator xi,
	//        BigInt Base pointer *B,
	//        The "carry" placeholder pointer
	//
	// Output: updates carry, and outputs the sign of number at position xi after adding carry.
	// Todo: handle the zero case--can't have a negative if the number is zero.
	// ----------------------------------------------------------------------------------------
	**xi += (long double)((long long)*carry);
	int sign = 1;
	if (**xi < 0) sign = -1;

	if (abs(**xi) > *B) {
		long double r;
		ddiv(xi,B,carry,&r);
		**xi = (long double)((long long)abs(r));
	} else *carry = 0;
	return sign;
}
/**
 * Implements Rearrange on a BigInt. Iterates 0 -> n.
 * @param X : the BigInt Pointer
 */
void rearrange(BigInt *X){
	long double carry = 0; int sign = 1; bool all_zeros = true;
	vector<long double>::iterator xi = X->number.begin();
	while(xi != X->number.end()){
		// carry carries through the sign. What matters is the sign of the last carry.
		sign *= rearrange_subproblem(&xi,&(X->base),&carry);
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
	return *signX * (**xiter) - *signY*(**yiter);
}
/**
 * Evaluates X + Y = Z, left to right, in an immutable fashion (creates new BigInt, Z)
 * @param X : Left hand BigInt Pointer
 * @param Y : Right hand BigInt Pointer
 * @return Z : solution to X + Y
 */
BigInt Add(BigInt *X,BigInt *Y){
	return SimpleArithmetic(X,Y,AddSubproblem);
}
/**
 * Evaluates X - Y = Z, left to right, in an immutable fashion (creates new BigInt, Z)
 * @param X : Left hand BigInt Pointer
 * @param Y : Right hand BigInt Pointer
 * @return Z : solution to X - Y, read left to right.
 */
BigInt Subtract(BigInt *X,BigInt *Y){
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

// ToDo: Test and document below
inline long double MultiplySubproblem(int *signX,vector<long double>::const_iterator *xiter,
                                      int *signY,vector<long double>::const_iterator *yiter){
	return (*signX * (**xiter)) * (*signY*(**yiter));
}
BigInt Multiply(BigInt *X,int q){
	BigInt Z(X->number.size(),X->base,1);

	vector<long double>::const_iterator xi = X->number.begin();
	vector<long double>::iterator zi = Z.number.begin();

	int abs_q = abs(q);
	int sign_q; if (q >= 0) sign_q = 1; else sign_q = -1;
	long double carry = 0;
	long double result;
	while(xi != X->number.end()){
		// solve first
		// ToDo: check if I need this trunc (if there are any rounding errors)
		result = (*xi+carry) * abs_q;
		// quotient is the carry, remainder is the value I keep in zi;
		ddiv(&result,&(X->base),&carry,&(*zi));
		zi++;xi++;
	}
	// adds a maximum of +1 to the length of the vector. can be postive or negative.
	// so there will be no "carry" on this operation.
	if (carry != 0){
		Z.number.push_back(carry);
	}
	// Finally, calculate the new sign of the BigInt:
	Z.sign *= sign_q;
	// All done, return Z.
	return Z;
}
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

// ToDo: Break up FFT code (written) complexity properly
template<typename long double>
vector<long double> FFT(vector<long double> *A,long double w){
	// Input: vector<long double> A length n, where long double is a number type, and w, where w is an nth root of unity in long double
	// Output: vector<long double> A, where values are A evaluated at n points.
	//
	// Base Case
	if (w == 1) return vector<long double> (1,1);
	// variables: int i, bool even
	ulong i; bool even = true; ulong n = A->size();
	// recursion variables:
	// ToDo: FFT clean up; pre-allocate, etc.
	vector<long double> ae; vector<long double> ao; vector<long double>::const_iterator Ai = A->begin();
	// Output variable
	vector<long double> O(n);

	// express A(x) in the form Ae(x^2) + x Ao(x2)
	// ae = {x_0, x_2, x_4, ..., x_n-2}
	// ao = {x_1, x_3, x_5, ..., x_n-1}
	while (Ai != A->end()){
		if (even){ ae.push_back(*Ai); even = false; }
		else     { ao.push_back(*Ai); even = true;  }
		Ai++;
	}

	// Divide
	vector<long double> Ae = FFT(&ae,w*w); // Evaluate A at even powers
	vector<long double> Ao = FFT(&ao,w*w); // Evaluate A at odd powers

	// Conquer (Merge)
	for (i = 0; i < n/2; i++){
		O[i] = Ae[i] + pow(w,i)*Ao[i];
		O[i + n/2] = Ae[i] - pow(w,i)*Ao[i];
	}
	return O;
}

// ToDo: FFT Multiply.
BigInt FFTMultiply(BigInt *A, BigInt *B){
	// Input: BigInt A, BigInt B
	// Output: A*B in O(nlgn)

	// ToDo: rearrange A, B in terms of highest minimum accurate base (or keep current base if it is sufficient!)


	// vars
	int i = 0;

	// 1) Pad out A and B to length >= max(A->size(),B->size()) >= 2 raised to some integer p
	//    ...small as possible!
	// ToDo: correct this to power of 2 >= 2*max(A->size(),B->size())
	long double n = (long double) max(A->number.size(),B->number.size());
	// padd out to n
	vector<long double> Ap(n,0);
	vector<long double> Bp(n,0);

	copy(A->number.begin(),A->number.end(),Ap.begin());
	copy(B->number.begin(),B->number.end(),Bp.begin());


	// next, calculate nth root of unity.
	long double Wn = -1*exp(2L * M_PIl / n);

	vector<long double> Av = FFT(&Ap,Wn);
	vector<long double> Bv = FFT(&Bp,Wn);

	// Multiply the values elementwise
	vector<long double> Zv(n);
	vector<long double>::const_iterator zvi = Zv.begin();
	vector<long double>::const_iterator avi = Av.begin();
	vector<long double>::const_iterator bvi = Bv.begin();
	while (avi != Av.end()){
		*zvi = *avi * *bvi;
		zvi++;avi++;bvi++;
	}
	// FFT result
	// Todo: get the sign of this thing figured out.
	BigInt Z(FFT(&Zv,-1L*Wn),A->base,1);
	// rearrange Z
	rearrange(&Z);

	// ToDo: if A,B accuracy rearrangement is done, then rearrange Z to have the old base!

	return Z;
}

// ToDo: SlowDivision (there must be a fast division.  This is hilariously absent in all the literature.)

int main() {
	// Tests: sign and accuracy. All big int vectors should be all positive.
	// ...sadly, until I can operate in ulonglong, this divides my total numbers territory.
	// ... but operating in ulonglong divides my speed, since this bars FFT.
	// ...but I need something that works.
	long double base = (long double) INT32_MAX;
	vector<long double> x_(10,(long double) base);BigInt x = BigInt(x_,base,1);
	//region 1) Add positive to positive - PASS
	vector<long double> y_(10,8.0);BigInt y = BigInt(y_,base,1);
    // BigInt z = Add(&x,&y);
	//endregion

	//region 2) Add positive to negative - PASS (also, vector is cleaned up).
	vector<long double> y2_(10,(long double)base);BigInt y2 = BigInt(y2_,base,-1);
    // BigInt z2 = Add(&y2,&x);


	// 3) Add negative to positive       - PASS
	BigInt z3 = Add(&x,&y2);


	// 4) Add negative to negative       - PASS
	vector<long double> x4(10,(long double) base); BigInt X4 = BigInt(x4,base,-1);
	vector<long double> y4(10,(long double) base); BigInt Y4 = BigInt(y4,base,-1);
    // BigInt Z4 = Add(&X4,&Y4);

	// 5) Add: test different size numbers PASS
	vector<long double> x5(5,(long double) base); BigInt X5 = BigInt(x5,base,1);
	vector<long double> y5(10,(long double) base);BigInt Y5 = BigInt(y5,base,1);
	BigInt Z5 = Add(&X5,&Y5);
	int xx = 1;
	// 5) Subtract small positive from large positive
	// todo: run test 5

	// 6) Subtract large positive from small positive

	// 7) Subtract large negative from small positive
	// 8) Subtract small negative from large positive

	// 9) Subtract large negative from small negative
	// 10)Subtract small negative from large negative

	std::cout << "Hello, World!" << std::endl;
	return 0;
}