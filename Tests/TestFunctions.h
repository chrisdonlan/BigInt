//
// Created by chris on 5/24/2017.
//

#ifndef BIGINT3_TESTFUNCTIONS_H
#define BIGINT3_TESTFUNCTIONS_H
#include <vector>
using namespace std;

vector<bool> epsilon_test(vector<long double> expected,vector<long double> observed,long double epsilon);
vector<vector<bool>> epsilon_test(vector<vector<long double>> expected,vector<vector<long double>> observed, long double epsilon);
void assert_true(vector<bool> test);
void assert_true(vector<vector<bool>> test);
#endif //BIGINT3_TESTFUNCTIONS_H
