//
// Created by Chris Donlan on 5/16/2017.
//

#include "BigIntTest.h"
#include <gmock/gmock.h>
#include "TestFunctions.h"
static long double epsilon = 1e-10;
using testing::Eq;
namespace {
	class ClassDeclaration: public testing::Test {
	public:
//		BigInt object;
		ClassDeclaration(){
//			object; // put method here
		}
	};
	class BigIntLogic: public testing::Test{
	public:
//		BigInt object;
		BigIntLogic(){
//			object;
		}
	};
	class BigIntAddition:public testing::Test{
	public:
//		BigInt object;
		BigIntAddition(){
//			object;
		}
	};
	class BigIntSubtraction:public testing::Test{
	public:
//		BigInt object;
		BigIntSubtraction(){
//			object;
		}
	};
	class BigIntMultiplication:public testing::Test{
	public:
//		BigInt object;
		BigIntMultiplication(){
//			object;
		}
	};
	class BigIntDivision:public testing::Test{
	public:
		BigIntDivision(){};
	};
}

TEST_F(ClassDeclaration, default_test) {
	BigInt x = BigInt();
	bool t0 = x.base == INT32_MAX;
	bool t1 = x.number.size() == 0;
	bool t2 = x.sign == 1;

	ASSERT_EQ(t0,true);
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(ClassDeclaration,with_ulong){
	BigInt x = BigInt(5,5L,1);
	bool t0 = x.base == 5L;
	bool t1 = x.sign == 1;
	bool t2 = x.number.capacity() == 5;

	ASSERT_EQ(t0,true);
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(ClassDeclaration,with_vector) {
	vector<long double> num(1, 15);
	BigInt y = BigInt(num, 2, 1);
	// this should rebase the vector, and it should be
	// length 4 when done.
	// Further, all values should be 1s.

	ASSERT_EQ(y.base, 2);
	ASSERT_EQ(y.number.size(), 4);
	ASSERT_EQ(y.sign, 1);

	vector<long double>::const_iterator ni = y.number.begin();
	bool all_ones = true;
	while (ni != y.number.end()) {
		if (*ni != 1) all_ones = false;
		ni++;
	}
	ASSERT_EQ(all_ones, true);
}
TEST_F(ClassDeclaration,with_vector2){
	vector<long double> bin_test{1L, 1L, 0L, 0L, 1L};
	vector<long double> num2(1, 19);
	BigInt y2 = BigInt(num2,2,1);
	bool equal = true;
	vector<long double>::const_iterator i1 = bin_test.begin();
	vector<long double>::const_iterator i2 = y2.number.begin();

	if (bin_test.size() == y2.number.size()){
		while(i1 != bin_test.end()){
			if (*i1 != *i2) equal = false;
			i1++;i2++;
		}
	} else equal = false;
	ASSERT_EQ(equal,true);
}
TEST_F(ClassDeclaration,with_vector3){
	vector<long double> x(1,55L);

	BigInt X(x,10,1);
	vector<long double> bi_x = X.number;
	vector<long double> expected = {5,5};
	vector<bool> test = epsilon_test(expected,bi_x,epsilon);
	assert_true(test);
}

TEST_F(BigIntAddition,sX_plus_sY){
	// ToDo: enhance addition test
	vector<long double> x(2,5L);
	vector<long double> y(2,5L);

	BigInt X(x,10,1);
	BigInt Y(y,10,1);
	BigInt Z = X + Y;

	// Base 10, answer should be: 110 reversed...
	long double a_[] = {0L,1L,1L};
	vector<long double> a = {0L,1L,1L};
	bool right_sign = 1 == Z.sign;
	bool right_value = Z.equals(a);
	ASSERT_EQ(right_value,true);
	ASSERT_EQ(right_sign,true);
}
TEST_F(BigIntSubtraction,sX_sub_sY){
	// ToDo: enhance subtraction test
	vector<long double> x(2,5L);
	vector<long double> y(2,5L);

	BigInt X(x,10,-1);
	BigInt Y(y,10,1);
	BigInt Z = X - Y;
	// Base 10, answer should be: 0 1 1 , sign = -1;
	long double a_[] = {0L,1L,1L};
	vector<long double> a = {0L,1L,1L};
	bool right_sign = -1 == Z.sign;
	bool right_value = Z.equals(a);
	ASSERT_EQ(Z.equals(a),true);
	ASSERT_EQ(right_sign,true);
}

TEST_F(BigIntLogic,X_equals_Y){
	// ToDo: Testing:X_equals_Y: make this test more robust.
	vector<long double> x(2,5L);
	vector<long double> y(2,5L);
	vector<long double> z(3,5L);
	vector<long double> alpha(1,5L);


	BigInt X(x,10,1);
	BigInt Y(y,10,1);
	BigInt Z(z,10,1);
	BigInt Alpha(alpha,10,1);
	BigInt Beta;
	BigInt Delta;

    bool equality1 = X.equals(&Y); // true
	bool equality2 = X.equals(&Z); // false
	bool equality3 = X.equals(&Alpha); // false
	bool equality4 = X.equals(&Beta);  // false
	bool equality5 = Beta.equals(&Delta); // true


	ASSERT_EQ(equality1,true);
	ASSERT_EQ(equality2,false);
	ASSERT_EQ(equality3,false);
	ASSERT_EQ(equality4,false);
	ASSERT_EQ(equality5,true);
}

TEST_F(BigIntLogic,X_gt_Y_1){
	vector<long double> x(2,5L);

	// X less than y through sign, otherwise equal
	BigInt X(x,10,-1);
	BigInt Y(x,10,1);

	bool gt = X > Y;
	bool lt = X < Y;
	ASSERT_EQ(gt,false);
	ASSERT_EQ(lt,true);
}
TEST_F(BigIntLogic,X_gt_Y_2){
	vector<long double> x(2,5L);
	vector<long double> y(3,5L);
	// This time, X is less than Y by virtue of length
	BigInt pX(x,10,1);
	BigInt pY(y,10,1);

	BigInt nX(x,10,-1);
	BigInt nY(y,10,-1);

	bool ppgt = pX > pY; ASSERT_EQ(ppgt,false);//false
	bool pngt = pX > nY; ASSERT_EQ(pngt,true);//true
	bool pplt = pX < pY; ASSERT_EQ(pplt,true);//true
	bool pnlt = pX < nY; ASSERT_EQ(pnlt,false);//false

	bool nngt = nX > nY; ASSERT_EQ(nngt,true);//true
	bool npgt = nX > pY; ASSERT_EQ(npgt,false);//false
	bool nnlt = nX < nY; ASSERT_EQ(nnlt,false);//false
	bool nplt = nX < pY; ASSERT_EQ(nplt,true);//true
}
TEST_F(BigIntLogic,X_gele_Y){
	//ToDo:BigIntLogic:X_ge_Y: improve this test
	vector<long double> x(2,5L);
	vector<long double> y(1,55L);

	BigInt pX(x,10,1);
	BigInt pY(y,10,1);

	BigInt nX(x,10,-1);
	BigInt nY(y,10,-1);

	bool ppgt = pX >= pY; ASSERT_EQ(ppgt,true);//true
	bool pngt = pX >= nY; ASSERT_EQ(pngt,true);//true
	bool pplt = pX <= pY; ASSERT_EQ(pplt,true);//true
	bool pnlt = pX <= nY; ASSERT_EQ(pnlt,false);//false

	bool nngt = nX >= nY; ASSERT_EQ(nngt,true);//true
	bool npgt = nX >= pY; ASSERT_EQ(npgt,false);//false
	bool nnlt = nX <= nY; ASSERT_EQ(nnlt,true);//true
	bool nplt = nX <= pY; ASSERT_EQ(nplt,true);//true
}


TEST_F(BigIntMultiplication,pn_bigint_times_p_num){
	vector<long double> x(2,5L);
	BigInt X(x,10,1);
	BigInt X2(x,10,-1);

	int c = 2;
	BigInt Y = X*c;
	BigInt Y2 = X2*c;

	vector<long double> r = {0,1,1};
	bool t1 = Y.equals(r);
	bool ysign = Y.sign == 1;

	bool t2 = Y2.equals(r);
	bool y2sign = Y2.sign == -1;
	ASSERT_EQ(t1 && ysign && t2 && y2sign,true);
}
TEST_F(BigIntMultiplication,pn_bigint_times_n_num){
	vector<long double> x(2,5L);
	BigInt X(x,10,1);
	BigInt X2(x,10,-1);

	int c = -2;
	BigInt Y = X*c;
	BigInt Y2 = X2*c;

	vector<long double> r = {0,1,1};
	bool t1 = Y.equals(r);
	bool t2 = Y2.equals(r);

	bool ysign = Y.sign == -1;
	bool y2sign = Y2.sign == 1;

	ASSERT_EQ(t1&&t2&&ysign&&y2sign,true);
}
TEST_F(BigIntMultiplication,pn_bigint_times_zero){
	vector<long double> x(2,5L);
	BigInt X(x,10,1);
	BigInt X2(x,10,-1);

	int c = 0;
	BigInt Y = X*c;
	BigInt Y2 = X2*c;

	bool t1 = Y.number.size() == 0;
	bool t2 = Y2.number.size() == 0;

	bool ys = Y.sign == 1;
	bool y2s = Y2.sign == 1;

	ASSERT_EQ(t1&&t2&&ys&&y2s,true);
}
TEST_F(BigIntMultiplication,pnbigint_times_pnbigint){
	vector<long double> x(1,11L);
	vector<long double> y(1,11L);
	BigInt X (x,10, 1);
	BigInt X2(x,10,-1);
	BigInt Y (y,10, 1);
	BigInt Y2(y,10,-1);


	BigInt Z = X*Y; // positive
	BigInt Z2 = X*Y2; // negative
	BigInt Z3 = X2*Y; // negative
	BigInt Z4 = X2*Y2; // positive

	vector<long double> t = {1L,2L,1L};

	bool t1 = Z.equals(t);
	bool t2 = Z2.equals(t);
	bool t3 = Z3.equals(t);
	bool t4 = Z4.equals(t);

	bool s1 = Z.sign ==  1;
	bool s2 = Z2.sign == -1;
	bool s3 = Z3.sign == -1;
	bool s4 = Z4.sign ==  1;

	bool test = t1&&t2&&t3&&t4&&s1&&s2&&s3&&s4;
	ASSERT_EQ(test,true);
}
TEST_F(BigIntMultiplication,pnbigint_times_pnbigint_diff_bases){
	vector<long double> x(1,11L);
	vector<long double> y(1,11L);
	BigInt X (x,10, 1);
	BigInt X2(x,10,-1);
	BigInt Y (y,12, 1);
	BigInt Y2(y,12,-1);

	// fixme: segmentation fault
	BigInt Z = X*Y;
	BigInt Z2 = X*Y2;
	BigInt Z3 = X2*Y;
	BigInt Z4 = X2*Y2;

	vector<long double> t = {1L,2L,1L};

	bool t1 = Z.equals(t);
	bool t2 = Z2.equals(t);
	bool t3 = Z3.equals(t);
	bool t4 = Z4.equals(t);

	bool s1 = Z.sign ==  1;
	bool s2 = Z2.sign == -1;
	bool s3 = Z3.sign == -1;
	bool s4 = Z4.sign ==  1;

	bool test = t1&&t2&&t3&&t4&&s1&&s2&&s3&&s4;
	ASSERT_EQ(test,true);
}
TEST_F(BigIntMultiplication,bigint_times_empty_bigint){
	vector<long double> x;
	vector<long double> y(1,11L);

	BigInt X(x,10,1);
	BigInt X2(x,10,-1);
	BigInt Y(y,10,1);

	BigInt Z = X*Y;

	bool t1 = Z.number.size() == 0;
	bool t2 = Z.sign == 1;

	ASSERT_EQ(t1&&t2,true);
}

TEST_F(BigIntDivision,pn_bigint_divided_by_num){
	vector<long double> x(2,6L);
	BigInt X(x,10,1);
	BigInt Y(x,10,-1);

	int c = 2;
	BigInt Z1 = X / c;
	BigInt Z2 = Y / c;

	vector<long double> r = {3,3};
	bool t1 = Z1.equals(r);
	bool t2 = Z1.sign == 1;

	bool t3 = Z2.equals(r);
	bool t4 = Z2.sign == -1;

	ASSERT_EQ(t1&&t2&&t3&&t4,true);
}