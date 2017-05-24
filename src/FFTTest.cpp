//
// Created by Chris Donlan on 5/23/2017.
//

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

