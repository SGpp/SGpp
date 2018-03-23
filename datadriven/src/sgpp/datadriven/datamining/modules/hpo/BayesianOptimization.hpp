/*
 * BayesianOptimization.hpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include "BOConfig.hpp"


namespace sgpp {
namespace datadriven {

class BayesianOptimization {
public:
  explicit BayesianOptimization(const std::vector<BOConfig>& initialConfigs);

  double acquisitionOuter(const base::DataVector& inp);

	//double mean(base::DataVector knew);
	//double var(base::DataVector knew, double kself);
	void updateGP();
	//double acquisitionPI(base::DataVector knew, double kself, double bestsofar);
	double acquisitionEI(base::DataVector knew, double kself, double bestsofar);
  void CholeskyDecomposition();
  void solveCholeskySystem(base::DataVector& x);
  BOConfig* main(BOConfig& prototype);
  double transformScore(double original);
  double kernel(double distance);



    // double expkernel(base::DataVector x1, base::DataVector x2);
protected:
	base::DataMatrix kernelmatrix;
  base::DataMatrix gleft;
  base::DataVector transformedOutput;
  base::DataVector rawScores;
  double bestsofar;
	bool screwedvar;
	double maxofmax;


  std::vector<BOConfig> allConfigs;
  // base::DataMatrix gright;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_ */
