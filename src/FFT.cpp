//
// Created by Chris Donlan on 5/22/2017.
//
#include <iostream>
#include "FFT.h"
vector<long double> nth_root_of_unity(long double n){
	return {1,2*M_PIl / n};
}
vector<long double> inverse_nth_root_of_unity(long double n){
	return {1,-2*M_PIl/n};
}
vector<vector<long double>> FFTPrep(vector<vector<long double>> *A){
	long n = 1;
	while(n < A->size()) n*=2;
	vector<vector<long double>> B(n,{0,0});
#pragma omp parallel for
	for(int i = 0; i < A->size(); i++){
		B[i][0] = (*A)[i][0];B[i][1] = (*A)[i][1];
	}
	return B;
}
vector<vector<long double>> FFTPrep(vector<long double> *A){
	long n = 1;
	while(n < A->size()) n*=2;
	vector<vector<long double>> B(n,{0,0});
#pragma omp parallel for
	for(int i = 0; i < A->size(); i++){
		B[i][0] = (*A)[i];
	}
	return B;
}
void PadTo(vector<vector<long double>> *A,long double target_size){
	while(A->size() < target_size) A->push_back({0,0});
}
vector<vector<long double>> PadCopy(vector<vector<long double>> *A,long double target_size){
	try{if(target_size < A->size()) throw 0;}catch (int i){cout << "target size must be greater than or equal to vector size.";}
	vector<vector<long double>> B(target_size,{0,0});
#pragma omp parallel for
	for (int i=0; i < A->size(); i++){
		B[i][0] = (*A)[i][0];B[i][1] = (*A)[i][1];
	}
	return B;
}
void PadTo(vector<long double> *A,long double target_size){
	while(A->size() < target_size) A->push_back(0);
}
vector<long double> PadCopy(vector<long double> *A,long double target_size){
	try{if(target_size < A->size()) throw 0;}catch (int i){cout << "target size must be greater than or equal to vector size.";}
	vector<long double> B(target_size,0);
#pragma omp parallel for
	for (int i = 0; i < A->size(); i++){
		B[i] = (*A)[i];
	}
}
void FFTMultiPrep(vector<long double> *A,vector<long double>*B){
	// Null case:
	if (A->size() == 0 && B->size() == 0){
		A->push_back(0);
		B->push_back(0);
	}
	while(A->size() != B->size()){
		if (A->size() < B->size()) A->push_back(0);
		else                       B->push_back(0);
	}

	long n = 1; long minimum = A->size()*2;
	while (n < minimum) n*=2;
	PadTo(A,n);
	PadTo(B,n);
}
void FFTMultiPrep(vector<vector<long double>> *A,vector<vector<long double>> *B) {
	// Null case:
	if (A->size() == 0 && B->size() == 0) {
		A->push_back({0,0});
		B->push_back({0,0});
	}

	// first, sizes have to be equal
	while(A->size() != B->size()){
		if (A->size() < B->size()) A->push_back({0,0});
		else                       B->push_back({0,0});
	}

	// next, calculate target
	long n = 1;long minimum = A->size()*2;
	while(n < minimum) n*=2;
	PadTo(A,n);
	PadTo(B,n);
}
vector<vector<long double>> FFT(vector<vector<long double>> *A,vector<long double>*w){
	/*
	 * Inputs: Vector of complex values A, B; polar form; complex values are vectors length 2: R, C
	 *
	 *         - assumed to be in standard order:  A0 + A1*x + A2*x^2 + ...
	 *
	 * Output: the transformed vector, Av (or Ap), Polar form
	 */
	// Step 1: return base case
	vector<long double> cw = cartesian(w);
	if (cw[0]==1) return *A;

	// Step 2: initialize variables
	long i; bool even = true; long n = A->size();
	vector<vector<long double>>::iterator Ai = A->begin();
	vector<vector<long double>> Out(n);
	vector<vector<long double>> ae, ao, Ae, Ao;

	// Algorithm: Divide and Conquer
	// Divide Part 1: Divide the polynomial between even and odd
	while(Ai != A->end()){
		if (even){ae.push_back(*Ai); even = false;}
		else     {ao.push_back(*Ai); even = true;}
		Ai++;}

	// Divide Part 2: Recurse
	//       calculate w^2 (squared root of unity)
	vector<long double> wsquared = {1,2*(*w)[1]};
	Ae = FFT(&ae,&wsquared);
	Ao = FFT(&ao,&wsquared);


	// Conquer
//#pragma omp parallel for
	for(i = 0; i < n/2; i++){
		vector<long double> wi      = {(*w)[0],(*w)[1]*i};
		vector<long double> xAo     = multiply_polar(&wi,&Ao[i]);
		vector<long double> cAe     = cartesian(&Ae[i]);
		vector<long double> cxAo    = cartesian(&xAo);
		vector<long double> cfirst  = add_cartesian(&cAe,&cxAo);
		vector<long double> csecond = subtract_cartesian(&cAe,&cxAo);


		Out[i]     = polar(&cfirst);
		Out[i+n/2] = polar(&csecond);
	}

	return Out;
}
vector<vector<long double>> FFTMultiplyComplex(vector<vector<long double>> A,vector<vector<long double>> B){
	/*
	 * Inputs: 2 vectored numbers in complex form: A0,A1,...,An;   B0,B1,...,Bm;
	 *
	 * Outputs: A*B, in complex form
	 *
	 * Assumptions:
	 *
	 *      - Complex Numbers: Assumed to be in POLAR COORDINATES!!!
	 *
	 *      - Length: vectors are *NOT* assumed to be equal length, a power of 2, 2n
	 *        **** In other words: vectors are assumed to be unmodified.
	 *
	 *      - Order: vectors are assumed to be a single number, of a given base B such that
	 *               A[0] = c0, A[1] = c1, ..., A[n] = cn in the polynomial:
	 *                    c0x^0 + c1x + c2x^2 + ... + cnx^n = A(x)
	 *
	 *               if, given traditional number, 1234 read left to right, in terms of the polynomial
	 *               4*10^0 + 3*10^1 + 2*10^2 + 1*10^3, the vector is assumed to be the reverse.
	 *
	 *      - Base B: Base B is assumed to be sufficiently small to guarantee the accuracy of the FFT
	 *                multiplication, which has error proportional to alpha = O(n*epsilon*B^2),
	 *                and alpha <= 6*n^2*B^2*log(n)*epsilon, where n is length of the vector,
	 *                B is the bit size?
	 */
	// Step 1: prep vectors--pad out to 2^p, such that 2^p >= 2*length of the longest vector
	FFTMultiPrep(&A,&B); long double n = (long double) A.size();
	vector<long double> w    = nth_root_of_unity(n);
	vector<long double> winv = inverse_nth_root_of_unity(n);

	// Step 4: FFT of vectors
	vector<vector<long double>> Av = FFT(&A,&w);
	vector<vector<long double>> Bv = FFT(&B,&w);

	// Elementwise Multiplication
	vector<vector<long double>> Cv = multiply_polar(&Av, &Bv);

	// Step 5: Get polynomial coefficients Todo: multiple of O(n) optimization possible here
	vector<vector<long double>> Cpi = FFT(&Cv,&winv);
	long double correction_factor = 1.0/n;
	vector<vector<long double>> Cp = multiply_polar(&correction_factor,&Cpi);
	return Cp;
}
vector<vector<long double>> FFTMultiplyComplex(vector<vector<long double>> *A,vector<vector<long double>> *B){
	/*
	 * Inputs: 2 vectored numbers in complex form: A0,A1,...,An;   B0,B1,...,Bm;
	 *
	 * Outputs: A*B, in complex form
	 *
	 * Assumptions:
	 *
	 *      - Complex Numbers: Assumed to be in POLAR COORDINATES!!!
	 *
	 *      - Length: vectors are *NOT* assumed to be equal length, a power of 2, 2n
	 *        **** In other words: vectors are assumed to be unmodified.
	 *
	 *      - Order: vectors are assumed to be a single number, of a given base B such that
	 *               A[0] = c0, A[1] = c1, ..., A[n] = cn in the polynomial:
	 *                    c0x^0 + c1x + c2x^2 + ... + cnx^n = A(x)
	 *
	 *               if, given traditional number, 1234 read left to right, in terms of the polynomial
	 *               4*10^0 + 3*10^1 + 2*10^2 + 1*10^3, the vector is assumed to be the reverse.
	 *
	 *      - Base B: Base B is assumed to be sufficiently small to guarantee the accuracy of the FFT
	 *                multiplication, which has error proportional to alpha = O(n*epsilon*B^2),
	 *                and alpha <= 6*n^2*B^2*log(n)*epsilon, where n is length of the vector,
	 *                B is the bit size?
	 */
	// Step 1: prep vectors--pad out to 2^p, such that 2^p >= 2*length of the longest vector
	FFTMultiPrep(A,B); long double n = (long double) A->size();
	vector<long double> w    = nth_root_of_unity(n);
	vector<long double> winv = inverse_nth_root_of_unity(n);

	// Step 4: FFT of vectors
	vector<vector<long double>> Av = FFT(A,&w);
	vector<vector<long double>> Bv = FFT(B,&w);

	// Elementwise Multiplication
	vector<vector<long double>> Cv = multiply_polar(&Av, &Bv);

	// Step 5: Get polynomial coefficients Todo: multiple of O(n) optimization possible here
	vector<vector<long double>> Cpi = FFT(&Cv,&winv);
	long double correction_factor = 1.0/n;
	vector<vector<long double>> Cp = multiply_polar(&correction_factor,&Cpi);
	return Cp;
}
vector<long double> FFTMultiply(vector<long double> *A,vector<long double> *B){
	vector<vector<long double>> Ac = Real2Complex(A);
	vector<vector<long double>> Bc = Real2Complex(B);
	vector<vector<long double>> Cc = FFTMultiplyComplex(&Ac,&Bc);
	return Complex2RealPolar(&Cc);
}
long double optimal_base(long size,long double base){
	long double precision = 1e-18;
	long double n = (long double) size;
	long double accuracy = base*base*precision*size;
	do {
		if (accuracy < 0.1){
			base*=2; n = fmaxl(1,n/ 2);
			accuracy = base*base*precision*n;
		} else {
			base /=2; n *= 2;
			accuracy = base*base*precision*n;
		}
	} while(accuracy < 0.1 || accuracy > 0.25);
	return base;
}
