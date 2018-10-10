// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp>

class Rand_double {
 public:
  Rand_double(double low, double high)
      : r(std::bind(std::uniform_real_distribution<>(low, high), std::default_random_engine())) {}

  double operator()() { return r(); }

 private:
  std::function<double()> r;
};

namespace sgpp {
namespace optimization {

void ASMatrixGradientMC::createMatrixMonteCarlo(size_t numPoints) {
  //  RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  Rand_double rd{0, 1};
  for (size_t i = 0; i < numPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    // todo (rehmemk) somehow the randomnumbergenerator results are much worse than the Rand_double
    // result
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    //    for (size_t d = 0; d < numDim; ++d) {
    //      randomVector[d] = rd();
    //    }
    //    std::cout << randomVector.toString() << std::endl;
    sgpp::base::DataVector gradient(numDim);
    objectiveFuncGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numPoints);
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
