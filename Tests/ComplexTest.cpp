//
// Created by Chris Donlan on 5/23/2017.
//

#include "TestFunctions.h"
#include "ComplexTest.h"
#include "../src/Complex.h"
#include <gmock/gmock.h>
#include <vector>


using namespace std;
using testing::Eq;
namespace{
	class PolarFunctions: public testing::Test {
	public:
		PolarFunctions(){}
	};
	class PolarVectorFunctions: public testing::Test {
	public:
		PolarVectorFunctions(){}
	};
	class CartesianFunctions: public testing::Test {
	public:
		CartesianFunctions(){}
	};
	class CartesianVectorFunctions: public testing::Test{
	public:
		CartesianVectorFunctions(){}
	};
}
long double truncto(long double *v,int decimal_place){
	long double raiseto = powl(10,(long double) decimal_place);
	*v *= raiseto;
	*v = truncl(*v);
	*v /= raiseto;
}
TEST_F(PolarFunctions,cartesian2polar_q1){
	vector<long double> c = {sqrtl(3)/2,0.5};
	vector<long double> c2= {0.5,sqrtl(3)/2};

	vector<long double> p = polar(&c);
	vector<long double> p2= polar(&c2);
	long double theta = M_PIl / 6;
	long double theta2 = 2*M_PIl/6;
	long double r = 1;
	long double epsilon = 1e-10;

	bool t1 = fabsl(r-p[0])<epsilon;
	bool t2 = fabsl(theta - p[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r-p2[0])<epsilon;
	t2 = fabsl(theta2-p2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(PolarFunctions,cartesian2polar_q2){
	vector<long double> c = {-sqrtl(3)/2,0.5};
	vector<long double> c2= {-0.5,sqrtl(3)/2};

	vector<long double> p = polar(&c);
	vector<long double> p2= polar(&c2);
	long double theta = 5*M_PIl / 6;
	long double theta2 = 4*M_PIl/6;
	long double r = 1;
	long double epsilon = 1e-10;

	bool t1 = fabsl(r-p[0])<epsilon;
	bool t2 = fabsl(theta - p[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r-p2[0])<epsilon;
	t2 = fabsl(theta2-p2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(PolarFunctions,cartesian2polar_q3){
	vector<long double> c = {-sqrtl(3)/2,-0.5};
	vector<long double> c2= {-0.5,-sqrtl(3)/2};

	vector<long double> p = polar(&c);
	vector<long double> p2= polar(&c2);
	long double theta = 7*M_PIl / 6;
	long double theta2 = 8*M_PIl/6;
	long double r = 1;
	long double epsilon = 1e-10;

	bool t1 = fabsl(r-p[0])<epsilon;
	bool t2 = fabsl(theta - p[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r-p2[0])<epsilon;
	t2 = fabsl(theta2-p2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(PolarFunctions,cartesian2polar_q4){
	vector<long double> c = {sqrtl(3)/2,-0.5};
	vector<long double> c2= {0.5,-sqrtl(3)/2};

	vector<long double> p = polar(&c);
	vector<long double> p2= polar(&c2);
	long double theta = -1*M_PIl / 6;
	long double theta2 = -2*M_PIl/6;
	long double r = 1;
	long double epsilon = 1e-10;

	bool t1 = fabsl(r-p[0])<epsilon;
	bool t2 = fabsl(theta - p[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r-p2[0])<epsilon;
	t2 = fabsl(theta2-p2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}

TEST_F(CartesianFunctions,polar2cartesian_q1){
	vector<long double> p  = {1,M_PIl / 6};
	vector<long double> p2 = {1,2*M_PIl / 6};
	vector<long double> c = cartesian(&p);
	vector<long double> c2 = cartesian(&p2);

	long double r = sqrtl(3)/2;
	long double i = 0.5;

	long double r2 = 0.5;
	long double i2 = sqrtl(3)/2;

	long double epsilon = 1e-10;

	bool t1,t2;
	t1 = fabsl(r-c[0])<epsilon;
	t2 = fabsl(i-c[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r2 - c2[0])<epsilon;
	t2 = fabsl(i2 - c2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(CartesianFunctions,polar2cartesian_q2){
	vector<long double> p  = {1,5*M_PIl / 6};
	vector<long double> p2 = {1,4*M_PIl / 6};
	vector<long double> c = cartesian(&p);
	vector<long double> c2 = cartesian(&p2);

	long double r = -sqrtl(3)/2;
	long double i = 0.5;

	long double r2 = -0.5;
	long double i2 = sqrtl(3)/2;

	long double epsilon = 1e-10;

	bool t1,t2;
	t1 = fabsl(r-c[0])<epsilon;
	t2 = fabsl(i-c[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r2 - c2[0])<epsilon;
	t2 = fabsl(i2 - c2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(CartesianFunctions,polar2cartesian_q3){
	vector<long double> p  = {1,7*M_PIl / 6};
	vector<long double> p2 = {1,8*M_PIl / 6};
	vector<long double> c = cartesian(&p);
	vector<long double> c2 = cartesian(&p2);

	long double r = -sqrtl(3)/2;
	long double i = -0.5;

	long double r2 = -0.5;
	long double i2 = -sqrtl(3)/2;

	long double epsilon = 1e-10;

	bool t1,t2;
	t1 = fabsl(r-c[0])<epsilon;
	t2 = fabsl(i-c[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r2 - c2[0])<epsilon;
	t2 = fabsl(i2 - c2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(CartesianFunctions,polar2cartesian_q4){
	vector<long double> p  = {1,11*M_PIl / 6};
	vector<long double> p2 = {1,10*M_PIl / 6};
	vector<long double> c = cartesian(&p);
	vector<long double> c2 = cartesian(&p2);

	long double r = sqrtl(3)/2;
	long double i = -0.5;

	long double r2 = 0.5;
	long double i2 = -sqrtl(3)/2;

	long double epsilon = 1e-10;

	bool t1,t2;
	t1 = fabsl(r-c[0])<epsilon;
	t2 = fabsl(i-c[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);

	t1 = fabsl(r2 - c2[0])<epsilon;
	t2 = fabsl(i2 - c2[1])<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}


TEST_F(CartesianFunctions,back_and_forth_c2p2c){
	bool t1,t2; long double epsilon = 1e-10;
	vector<long double> c = {1,2};

	vector<long double> cp = polar(&c);
	vector<long double> cc = cartesian(&cp);

	t1 = fabsl(cc[0]-c[0])<epsilon;
	t2 = fabsl(cc[1]-c[1])<epsilon;

	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(CartesianFunctions,subtract_cartesian){
	vector<long double> c1 = {1,2};
	vector<long double> c2 = {0.5,1};
	vector<long double> c3 = subtract_cartesian(&c1,&c2);

	ASSERT_EQ(c3[0] == 0.5 && c3[1] == 1,true);
}
TEST_F(CartesianFunctions,add_cartesian){
	vector<long double> c1 = {1,2};
	vector<long double> c2 = {0.5,1};
	vector<long double> c3 = add_cartesian(&c1,&c2);

	ASSERT_EQ(c3[0] == 1.5 && c3[1] == 3,true);
}

TEST_F(PolarFunctions,back_and_forth_p2c2p){
	bool t1,t2; long double epsilon = 1e-10;
	vector<long double> p = {1,M_PIl/6};

	vector<long double> pc = cartesian(&p);
	vector<long double> pp = polar(&pc);

	t1 = fabsl(pp[0]-p[0])<epsilon;
	t2 = fabsl(pp[1]-p[1])<epsilon;

	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(PolarFunctions,subtract_polar){
	bool t1,t2;
	vector<long double> c1 = {1,2};
	vector<long double> c2 = {0.5,1};

	vector<long double> p1 = polar(&c1);
	vector<long double> ca = cartesian(&p1);
	vector<long double> p2 = polar(&c2);
	vector<long double> pc = subtract_polar(&p1,&p2);

	vector<long double> cc = cartesian(&pc);

	long double epsilon = 1e-10;

	t1 = fabsl(cc[0]-0.5)<epsilon;
	t2 = fabsl(cc[1]-1)<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}
TEST_F(PolarFunctions,add_polar){
	bool t1,t2; long double epsilon = 1e-10;
	vector<long double> c1 = {1,2};
	vector<long double> c2 = {0.5,1};

	vector<long double> p1 = polar(&c1);
	vector<long double> p2 = polar(&c2);
	vector<long double> p3 = add_polar(&p1,&p2);

	vector<long double> cc = cartesian(&p3);

	t1 = fabsl(cc[0]-1.5)<epsilon;
	t2 = fabsl(cc[1]-3)<epsilon;
	ASSERT_EQ(t1,true);
	ASSERT_EQ(t2,true);
}


TEST_F(PolarFunctions,multiply_polars){
	long double epsilon = 1e-10;
	vector<long double> c1 = {1,2},c2 = {3,4},expected = {-5,10};
	vector<long double> p1 = polar(&c1),p2 = polar(&c2);

	vector<long double> pr = multiply_polar(&p1,&p2);
	vector<long double> result = cartesian(&pr);

	vector<bool> t = epsilon_test(expected,result,epsilon);
	for(int i = 0; i < t.size(); i++){
		ASSERT_EQ(t[i],true);
	}
}
TEST_F(CartesianFunctions,multiply_cartesians){
	long double epsilon = 1e-10;
	vector<long double> c1 = {1,2};
	vector<long double> c2 = {3,4};

	vector<long double> exp = {-5,10};

	vector<long double> cm = multiply_cartesian(&c1,&c2);

	vector<bool> t = epsilon_test(exp,cm,epsilon);
	for(int i = 0; i < t.size(); i++){
		ASSERT_EQ(t[i],true);
	}
}
TEST_F(PolarFunctions,scalar_times_polar){
	long double epsilon = 1e-10; bool t1,t2;
	vector<long double> p = {1,M_PIl/6};
	long double q = 4;

	vector<long double> result = multiply_polar(&q,&p);
	vector<long double> expected = {4,M_PIl/6};

	vector<bool> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(CartesianFunctions,scalar_times_cartesian){
	long double epsilon = 1e-10; bool t1,t2;
	vector<long double> c = {1,2};
	long double q = 4;

	vector<long double> result = multiply_cartesian(&q,&c);
	vector<long double> expected = {4,8};

	vector<bool> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}


// New Section: Vector Functions (acting on a vector of complex values)
TEST_F(PolarVectorFunctions,convert_cartesian_vector){
	long double epsilon = 1e-10;
	vector<vector<long double>> c1 = {{1,2},{3,4},{5,6}};
	vector<vector<long double>> expected = {polar(&(c1[0])),polar(&(c1[1])),polar(&(c1[2]))};

	vector<vector<long double>> p = polar(&c1);

	vector<vector<bool>> test = epsilon_test(expected,p,epsilon);
	assert_true(test);
}
TEST_F(CartesianVectorFunctions,convert_polar_vector){
	long double epsilon = 1e-10;
	vector<vector<long double>> expected = {{1,2},{3,4},{5,6}};
	vector<vector<long double>> p1 = {polar(&(expected[0])),polar(&(expected[1])),polar(&(expected[2]))};

	vector<vector<long double>> c = cartesian(&p1);

	vector<vector<bool>> test = epsilon_test(expected,c,epsilon);
	assert_true(test);
}
TEST_F(PolarVectorFunctions,hadamard_polar_product){
	long double epsilon = 1e-10;
	// elementwise multiplication
	vector<vector<long double>> c1 = {{1,2},{3,4},{5,6}};
	vector<vector<long double>> c2 = {{7,8},{9,10},{11,12}};

	vector<vector<long double>> p1 = polar(&c1),p2 = polar(&c2);
	vector<vector<long double>> expected(3);
	expected[0] = multiply_polar(&p1[0],&p2[0]);
	expected[1] = multiply_polar(&p1[1],&p2[1]);
	expected[2] = multiply_polar(&p1[2],&p2[2]);

	vector<vector<long double>> result = multiply_polar(&p1,&p2);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(CartesianVectorFunctions,hadamard_cartesian_product){
	long double epsilon = 1e-10;
	// elementwise multiplication
	vector<vector<long double>> c1 = {{1,2},{3,4},{5,6}};
	vector<vector<long double>> c2 = {{7,8},{9,10},{11,12}};

	vector<vector<long double>> expected(3);
	expected[0] = multiply_cartesian(&c1[0],&c2[0]);
	expected[1] = multiply_cartesian(&c1[1],&c2[1]);
	expected[2] = multiply_cartesian(&c1[2],&c2[2]);

	vector<vector<long double>> result = multiply_cartesian(&c1,&c2);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(PolarVectorFunctions,scalar_times_polar_vector){
	long double q = 2,epsilon = 1e-10;
	vector<vector<long double>> p = {{1,M_PIl/6},{2,M_PIl/8},{3,M_PIl/10}};
	vector<vector<long double>> expected = {{2,M_PIl/6},{4,M_PIl/8},{6,M_PIl/10}};

	vector<vector<long double>> result = multiply_polar(&q,&p);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(PolarVectorFunctions,scalar_vector_times_polar_vector){
	vector<long double> q = {2,3,4}; long double epsilon = 1e-10;
	vector<vector<long double>> p = {{1,M_PIl/6},{2,M_PIl/8},{3,M_PIl/10}};
	vector<vector<long double>> expected = {{2,M_PIl/6},{6,M_PIl/8},{12,M_PIl/10}};
	vector<vector<long double>> result = multiply_polar(&q,&p);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(CartesianVectorFunctions,scalar_times_cartesian_vector){
	long double q = 2,epsilon = 1e-10;
	vector<vector<long double>> p = {{1,2},{3,4},{5,6}};
	vector<vector<long double>> expected = {{2,4},{6,8},{10,12}};
	vector<vector<long double>> result = multiply_cartesian(&q,&p);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(CartesianVectorFunctions,scalar_vector_times_cartesian_vector){
	vector<long double> q = {2,3,4};long double epsilon = 1e-10;
	vector<vector<long double>> p = {{1,2},{3,4},{5,6}};
	vector<vector<long double>> expected = {{2,4},{9,12},{20,24}};
	vector<vector<long double>> result = multiply_cartesian(&q,&p);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}

TEST_F(PolarVectorFunctions,complex_to_real){
	long double epsilon = 1e-10;
	vector<vector<long double>> p = {{2,M_PIl/2},{3,M_PIl/3},{4,M_PIl/4}};
	vector<long double> expected(3);
	vector<long double> c1 = cartesian(&(p[0])); expected[0] = c1[0];
	vector<long double> c2 = cartesian(&(p[1])); expected[1] = c2[0];
	vector<long double> c3 = cartesian(&(p[2])); expected[2] = c3[0];

	vector<long double> result = Complex2RealPolar(&p);

	vector<bool> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(PolarVectorFunctions,real_to_complex){
	long double epsilon = 1e-10;
	vector<long double> real = {1,2,3};
	vector<vector<long double>> expected = {{1,0},{2,0},{3,0}};

	vector<vector<long double>> result = Real2Complex(&real);

	vector<vector<bool>> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}
TEST_F(CartesianVectorFunctions,complex_to_real){
	long double epsilon = 1e-10;
	vector<vector<long double>> p = {{2,3},{4,5},{6,7}};
	vector<long double> expected = {2,4,6};

	vector<long double> result = Complex2RealCartesian(&p);

	vector<bool> test = epsilon_test(expected,result,epsilon);
	assert_true(test);
}







