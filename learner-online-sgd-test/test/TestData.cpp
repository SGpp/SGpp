#include "TestData.hpp"

namespace sg {

	namespace test {

		TestData::TestData(int dim_, int level_, size_t trainSize_)
		{
			dim = dim_;
			level = level_;
			trainSize = trainSize_;

			trainData = new sg::base::DataMatrix (trainSize, (size_t)dim);
			classes = new sg::base::DataVector (trainSize);
		}

		sg::base::DataMatrix TestData::getTrainData()
		{
			return *trainData;
		}

		sg::base::DataVector TestData::getClasses()
		{
			return *classes;
		}

		TestData::~TestData()
		{
			if (trainData != NULL)
				delete trainData;

			if (classes != NULL)
				delete classes;
		}
	}
}
