#include "PerformanceTest.hpp"
#include <iostream>

namespace sg {

  namespace test {

    PerformanceTest::PerformanceTest(bool verbose_) : error(NULL), mse(0) {
	  verbose = verbose_;
	  sg::datadriven::LearnerRegularizationType r = sg::datadriven::Identity;
	  learner = new sg::datadriven::LearnerSGD (r, true, false);
	}

	void PerformanceTest::runFriedman(int dim, int level, size_t trainSize, int seed, size_t maxIterations, const double eps, const double lambda, const double gamma) {

		FriedmanTestData testData(dim, level, trainSize);
		testData.generate(seed);

		sg::base::RegularGridConfiguration gridConfig;
		gridConfig.dim_ = dim;
		gridConfig.level_ = level;
		gridConfig.type_ = sg::base::ModLinear;

		sg::base::DataMatrix trainData = testData.getTrainData();
		sg::base::DataVector classes = testData.getClasses();

		learner->train(trainData, classes, gridConfig, maxIterations, eps, lambda, gamma); 

		updateErrors(*(learner->getGrid()), *(learner->getAlpha()), trainData, classes);
	}

	void PerformanceTest::updateErrors(sg::base::Grid& grid, sg::base::DataVector& alpha, sg::base::DataMatrix& trainData, sg::base::DataVector& classes)
	{
		int trainSize = trainData.getNrows();

		// Error vector
		error = new sg::base::DataVector(trainSize);
		error->setAll(0.0);

		sg::base::DataVector result(trainSize);
		sg::op_factory::createOperationMultipleEval(grid, &trainData)->mult(alpha, result);	

		for (int i = 0; i < trainSize; i++)
			error->set(i, result[i] - classes[i]);

		// MSE
		double sum;
		for (int i = 0; i < trainSize; i++)
			sum += pow(fabs(error->get(i)), 2);

		mse = sum/trainSize;
	}

	sg::base::DataVector PerformanceTest::getError()
	{
		return sg::base::DataVector(*error);
	}

	double PerformanceTest::getMSE()
	{
		return mse;
	}

	PerformanceTest::~PerformanceTest()
	{
		if (error != NULL)
			delete error;

		if (learner != NULL)
			delete learner;
	}
  }
}
