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

	void updateGP();
	double acquisitionEI(base::DataVector knew, double kself, double bestsofar);
  void decomposeCholesky(base::DataMatrix &km, base::DataMatrix &gnew);
  void solveCholeskySystem(base::DataMatrix &gmatrix, base::DataVector &x);
  BOConfig* main(BOConfig& prototype);
  double transformScore(double original);
  double kernel(double distance);
  void fitScales();
  double likelihood(const base::DataVector& inp);




  // double expkernel(base::DataVector x1, base::DataVector x2);
protected:
	base::DataMatrix kernelmatrix;
  base::DataMatrix gleft;
  base::DataVector transformedOutput;
  base::DataVector rawScores;
  base::DataVector scales;
  double bestsofar;
	bool screwedvar;
	double maxofmax;


  std::vector<BOConfig> allConfigs;
  // base::DataMatrix gright;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_ */
