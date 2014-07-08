#include <random>
#include "FriedmanTestData.hpp"
#include "sgpp_base.hpp"

namespace sg {

	namespace test {
		
		FriedmanTestData::FriedmanTestData (int dim_, int level_, size_t trainSize_) : TestData (dim_, level_, trainSize_) {
		}

        void FriedmanTestData::generate(int seed)
		{
			double pi = 3.14159265359;

			std::default_random_engine generator;
			generator.seed(seed);
			std::uniform_real_distribution<double> d1(0.0, 100.0);
			std::uniform_real_distribution<double> d2(40 * pi, 560 * pi);
			std::uniform_real_distribution<double> d3(0.0, 1.0);
			std::uniform_real_distribution<double> d4(1.0, 11.0);
			std::normal_distribution<double> deps(0.0, 125.0);

			for (unsigned int i = 0; i < trainSize; i++) {
				double x1 = d1(generator);
				double x2 = d2(generator);
				double x3 = d3(generator);
				double x4 = d4(generator);
				double eps = deps(generator);
				double y = sqrt(x1 * x1 + (x2 * x3 - 1.0 / x2 * x4) * (x2 * x3 - 1.0 / x2 * x4)) + eps;
				trainData->set(i, 0, x1);
				trainData->set(i, 1, x2);
				trainData->set(i, 2, x3);
				trainData->set(i, 3, x4);
				classes->set(i, y);
			}

			for (int d = 0; d < dim; d++)
				trainData->normalizeDimension(d);
		}
	}
}
