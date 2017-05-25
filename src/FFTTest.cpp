//
// Created by Chris Donlan on 5/23/2017.
//
#include "../Tests/TestFunctions.h"
#include "FFTTest.h"
#include "FFT.h"
#include <gmock/gmock.h>
#include <vector>
using namespace std;
using testing::Eq;
namespace{
	class FFTTesting: public testing::Test{
	public:
		FFTTesting(){}
	};
	class FFTMultiplyComplexTesting: public testing::Test{
	public:
		FFTMultiplyComplexTesting(){}
	};
}
static long double epsilon = 1e-10;
TEST_F(FFTTesting,Test_FFT_trivial) {
	vector<vector<long double>> A = {{1,0}};
	vector<long double> nth_root = nth_root_of_unity(1);
	vector<vector<long double>> C = FFT(&A,&nth_root);

	vector<vector<long double>> expected = {{1,0}};
	vector<vector<bool>> test = epsilon_test(C,expected,epsilon);
	assert_true(test);
}
TEST_F(FFTTesting,Test_FFT_double) {
	vector<vector<long double>> A = {{1,0},{2,0}};
	vector<long double> nth_root = nth_root_of_unity(2);
	vector<long double> nth_root_inv = inverse_nth_root_of_unity(2);
	vector<vector<long double>> C = FFT(&A,&nth_root);
	vector<vector<long double>> D_ = FFT(&C,&nth_root_inv);
	long double q = 1./2.;
	vector<vector<long double>> D  = multiply_polar(&q,&D_);
	vector<vector<bool>> test = epsilon_test(D,A,epsilon);
	assert_true(test);
}
TEST_F(FFTTesting,Test_FFT_long) {
	vector<vector<long double>> A = {{1,0},{2,0},{3,0},{4,0}};
	vector<long double> nth_root = nth_root_of_unity(4);
	vector<long double> nth_root_inv = inverse_nth_root_of_unity(4);
	vector<vector<long double>> C = FFT(&A,&nth_root);
	vector<vector<long double>> D_ = FFT(&C,&nth_root_inv);
	long double q = 1./4.;
	vector<vector<long double>> D  = multiply_polar(&q,&D_);
	vector<vector<bool>> test = epsilon_test(D,A,epsilon);
	assert_true(test);
}
TEST_F(FFTTesting,Test_FFT_Setup){
	vector<vector<long double>> A = {{1,0},{2,0},{3,0},{4,0},{5,0}};
	vector<vector<long double>> A_ = FFTPrep(&A);
	// should be first power of 2 greater or equal length A_
	vector<vector<long double>> expected = {{1,0},{2,0},{3,0},{4,0},{5,0},{0,0},{0,0},{0,0}};

	vector<vector<bool>> etest = epsilon_test(A_,expected,epsilon);
	bool stest = A_.size() == expected.size();

	assert_true(etest);
	ASSERT_EQ(stest,true);
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_MultiplyComplex){
	// Consider (1 + 2x + 3x^2) * (4 + 5x + 6x^2)
	// = 4 + 5x + 6x^2 + 8x + 10x^2 + 12x^3 + 12x^2 + 15x^3 + 18x^4
	// = 4 + 5x + 8x + 6x^2 + 10x^2 + 12x^2 + 12x^3 + 15x^3 + 18x^4
	// = 4 + 13x + 28x^2 + 27x^3 + 18x^4

	vector<vector<long double>> A = {{1,0},{2,0},{3,0}};
	vector<vector<long double>> B = {{4,0},{5,0},{6,0}};
	vector<long double> expected = {4,13,28,27,18,0};
	vector<vector<long double>> C_ =FFTMultiplyComplex(A,B);
	vector<long double> C = Complex2RealPolar(&C_);

	vector<bool> test = epsilon_test(expected,C,epsilon);
	assert_true(test);
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multiply){
	// Consider (1 + 2x + 3x^2) * (4 + 5x + 6x^2)
	// = 4 + 5x + 6x^2 + 8x + 10x^2 + 12x^3 + 12x^2 + 15x^3 + 18x^4
	// = 4 + 5x + 8x + 6x^2 + 10x^2 + 12x^2 + 12x^3 + 15x^3 + 18x^4
	// = 4 + 13x + 28x^2 + 27x^3 + 18x^4

	vector<long double> A = {1,2,3};
	vector<long double> B = {4,5,6};
	vector<long double> expected = {4,13,28,27,18,0};
	vector<long double> C =FFTMultiply(&A,&B);
	vector<bool> test = epsilon_test(expected,C,epsilon);
	assert_true(test);
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multply_zeros_left){
	// Todo
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multiply_zeros_right){
	// Todo
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multiply_zeros){
	// Todo
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multiply_null_left){
	// Todo
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multiply_null_right){
	// Todo
}
TEST_F(FFTMultiplyComplexTesting,Test_FFT_Multiply_nulls){
	// Todo
}

// ToDo: test negatives
// ToDo: test randoms