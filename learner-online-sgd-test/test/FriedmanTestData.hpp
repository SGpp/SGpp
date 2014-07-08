#ifndef FRIEDMANTESTDATA_HPP
#define FRIEDMANTESTDATA_HPP

#include "TestData.hpp"

namespace sg {

	namespace test {

		class FriedmanTestData: public TestData {

			public:
				FriedmanTestData (int dim_, int level_, size_t trainSize_);
			    void generate (int seed);
		};
	}
}

#endif
