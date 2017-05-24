//
// Created by Chris Donlan on 5/22/2017.
//
#include "FFT.h"
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
	if (cw[0]==1) return {{1,0}};

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
#pragma omp parallel for
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
	// ToDo: finish FFTMultiplyComplex
	long double n=1,ni=max(A.size(),B.size()),p=0;

	// Step 1: calculate n to use: n is power of 2, and >= 2n
	while(n < 2*ni){n*=2;p++;}
	vector<long double> w    = {1,   2*M_PIl / n};
	vector<long double> winv = {1,-1*2*M_PIl / n};

	// Step 2: pad out to the proper length (power of 2 >= 2n)
	while(A.size() < n) A.push_back({0,0});
	while(B.size() < n) B.push_back({0,0});

	// Step 3: nth root of unity

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
