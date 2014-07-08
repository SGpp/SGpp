#ifndef PERFORMANCETEST_HPP
#define PERFORMANCETEST_HPP


#include "test/FriedmanTestData.hpp"

#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"

namespace sg {

	namespace test {

		class PerformanceTest {

			public:
				PerformanceTest (bool verbose_);

				void runFriedman (int dim, int level, size_t trainSize,int seed, size_t maxIterations, double eps, double lambda, double gamma);

				sg::base::DataVector getError();
				double getMSE();

				~PerformanceTest();

				// get trainData
				// get classes
				// etc.

			private:
				bool verbose;

				double mse;
				sg::base::DataVector* error;
				sg::datadriven::LearnerSGD* learner;

			 	void updateErrors(sg::base::Grid& grid, sg::base::DataVector& alpha, sg::base::DataMatrix& trainData, sg::base::DataVector& classes);
		};

	}
}

#endif /* LEARNERSGD_HPP */
