//
// Created by Chris Donlan on 5/16/2017.
//

#ifndef BIGINT3_BIGINT_H
#define BIGINT3_BIGINT_H
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class BigInt {
public:
	bool rearranged;
	vector<long double> number;
	long double base;
	int sign;
	BigInt(ulong size,long double B,int s);
	BigInt(vector<long double> num,long double B,int s);
	BigInt();
	bool equals(BigInt *X);
	bool equals(BigInt X);
	bool equals(vector<long double> *x);
	bool equals(vector<long double> x);
	void update();

	// overload ++ (both sides), += here
};
BigInt Add(BigInt *X,BigInt *Y);
BigInt Subtract(BigInt *X,BigInt *Y);
BigInt Multiply(BigInt *X,int q);
BigInt Multiply(BigInt *A, BigInt *B);
BigInt Divide(BigInt *X, int q,long double * remainder);

inline bool operator >(BigInt X,BigInt Y){
	X.update();Y.update();
	// ToDo:Optimize: BigInt > BigInt
	bool greater_than = false;
	// deal with both null size
	if (X.number.size() == 0 && Y.number.size() == 0) return false;

	// Pad up to same size
	while(X.number.size() != Y.number.size()){
		if (X.number.size() < Y.number.size()){
			X.number.push_back(0);
		} else Y.number.push_back(0);
	}

	vector<long double>::const_iterator xi = X.number.end();
	vector<long double>::const_iterator yi = Y.number.end();
	do{
		xi--;yi--;
		if (*xi != 0 || *yi !=0){
			if (*xi * X.sign == *yi * Y.sign){
				continue;
			}
			else if (*xi * X.sign > *yi * Y.sign)
			{
				greater_than = true;
				break;
			} else // *xi*X.sign < *yi * Y.sign
			{
				greater_than = false;
				break;
			}
		}
	}while (xi!= X.number.begin());

	// Clean up Padding
	xi = X.number.end();
	yi = Y.number.end();
	do {
		xi--;
		X.number.pop_back();
	} while(*xi == 0 && xi != X.number.begin());
	do {
		yi--;
		Y.number.pop_back();
	} while(*xi == 0 && xi != X.number.begin());

	return greater_than;
}
inline bool operator <(BigInt X,BigInt Y){
	return Y > X;
}
inline bool operator ==(BigInt X,BigInt Y){
	// ToDo:Optimize: BigInt == BigInt in the long winded form
	X.update();Y.update();
	return X.equals(Y);
}
inline bool operator <=(BigInt X,BigInt Y){
	// ToDo:Optimize: BigInt <= BigInt in the long winded form

	bool a = X< Y;
	bool b = X == Y;
	return X < Y || X == Y;
}
inline bool operator >=(BigInt X,BigInt Y){
	// ToDo:Optimize: BigInt >= BigInt in the long winded form
	return X > Y || X == Y;
}

// ToDo:Test +(BigInt X,BigInt Y)
inline BigInt operator +(BigInt X,BigInt Y){return Add(&X,&Y);}
// ToDo:Test -(BigInt X,BigInt Y)
inline BigInt operator -(BigInt X,BigInt Y){return Subtract(&X,&Y);}
// ToDo:Test *(BigInt X,BigInt Y)
inline BigInt operator *(BigInt X,BigInt Y){return Multiply(&X,&Y);}
// ToDo:Test *(BigInt X,int q)
inline BigInt operator *(BigInt X,int q){return Multiply(&X,q);}
// ToDo:Test *(int q,BigInt X)
inline BigInt operator *(int q,BigInt X){return Multiply(&X,q);}
// ToDo:Test /(BigInt X,int q)
inline BigInt operator /(BigInt X,int q){long double r; return Divide(&X,q,&r);}

// ToDo:Implement /(int q,BigInt X)
// ToDo:Test /(int q,BigInt X)
inline BigInt operator /(int q,BigInt X){try{throw 0;}catch (int e){cout << "q/BigInt not implemented yet.";}}

// ToDo:Implement /(int q,BigInt X)
// Todo:Test: /(BigInt X,BigInt Y)
inline BigInt operator /(BigInt X,BigInt Y){try{throw 0;}catch (int e){cout << "BigInt/BigInt not implemented yet.";}}



#endif //BIGINT3_BIGINT_H
