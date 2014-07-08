#ifndef TESTDATA_HPP
#define TESTDATA_HPP

#include "sgpp_base.hpp"

namespace sg {

	namespace test {

		class TestData {

			public:
				TestData(int dim_, int level_, size_t trainSize_);
				~TestData();
				sg::base::DataMatrix getTrainData();
				sg::base::DataVector getClasses();

			protected:
				int dim;
				int level;
				size_t trainSize;
				sg::base::DataMatrix* trainData;
				sg::base::DataVector* classes;
		};
	}
}

#endif
