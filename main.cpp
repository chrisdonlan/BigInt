#include <iostream>
#include <vector>
#include <cmath>
#include "src/BigInt.h"

// Google Test Includes
#include <gtest/gtest.h>
#include<gmock/gmock.h>
#include "src/BigInt.h"


using namespace std;

int main(int argc,char* argv[]) {
	// Tests: sign and accuracy. All big int vectors should be all positive.
	// ...sadly, until I can operate in ulonglong, this divides my total numbers territory.
	// ... but operating in ulonglong divides my speed, since this bars FFT.
	// ...but I need something that works.
//	long double base = (long double) INT32_MAX;
//	vector<long double> x_(10,(long double) base);BigInt x = BigInt(x_,base,1);
//	vector<long double> y_(10,8.0);BigInt Y = BigInt(y_,base,1);
//	vector<long double> y2_(10,(long double)base);BigInt y2 = BigInt(y2_,base,-1);
//	vector<long double> x4(10,(long double) base); BigInt X4 = BigInt(x4,base,-1);
//	vector<long double> y4(10,(long double) base); BigInt Y4 = BigInt(y4,base,-1);
//	vector<long double> x5(5,(long double) base); BigInt X5 = BigInt(x5,base,1);
//	vector<long double> y5(10,(long double) base);BigInt Y5 = BigInt(y5,base,1);

	//region 1) Add positive to positive - PASS

    // BigInt z = Add(&x,&Y);
	//endregion

	//region 2) Add positive to negative - PASS (also, vector is cleaned up).

//     BigInt z2 = Add(&y2,&x);
	// 3) Add negative to positive       - PASS
//	BigInt z3 = Add(&x,&y2);
	// 4) Add negative to negative       - PASS
//     BigInt Z4 = Add(&X4,&Y4);
	// 5) Add: test different size numbers PASS
//	BigInt Z5 = Add(&X5,&Y5);
	// 6) Subtract small positive from large positive
//	BigInt Z6 = Subtract(&X5,&Y);
//	int xx = 1;
	// 6) Subtract large positive from small positive

	// 7) Subtract large negative from small positive
	// 8) Subtract small negative from large positive

	// 9) Subtract large negative from small negative
	// 10)Subtract small negative from large negative

	testing::InitGoogleTest(&argc,argv);
	RUN_ALL_TESTS();
	std::cout << "Hello, World!" << std::endl;
	return 0;
}