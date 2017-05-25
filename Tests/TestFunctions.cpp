//
// Created by chris on 5/24/2017.
//

#include <gtest/gtest.h>
#include <math.h>
#include "TestFunctions.h"
vector<vector<bool>> epsilon_test(vector<vector<long double>> expected,vector<vector<long double>> observed, long double epsilon){
	int i;
	vector<vector<bool>> results;
	for(i = 0; i < expected.size(); i++){
		if(i >= observed.size()) results.push_back({false,false});
		else {
			bool a = fabsl(expected[i][0]-observed[i][0])<epsilon;
			bool b = fabsl(expected[i][1]-observed[i][1])<epsilon;
			results.push_back({a,b});
		}
	}
	return results;
}
vector<bool> epsilon_test(vector<long double> expected,vector<long double> observed,long double epsilon){
	int i;
	vector<bool> results;
	for(i = 0; i<expected.size();i++){
		if (i >= observed.size()) results.push_back(false);
		else results.push_back(fabsl(expected[i]-observed[i])<epsilon);
	}
	return results;
}

void assert_true(vector<vector<bool>> test){
	for(int i = 0; i < test.size(); i++){
		for (int j = 0; j < test[i].size(); j++){
			ASSERT_EQ(test[i][j],true);
		}
	}
}
void assert_true(vector<bool> test){
	for (int i= 0; i < test.size(); i++)
		ASSERT_EQ(test[i],true);
}