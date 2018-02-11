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



namespace sgpp {
namespace datadriven {

class BayesianOptimization {
public:
	BayesianOptimization(double firstvalue);
	double mean(base::DataVector knew);
	double var(base::DataVector knew, double kself);
	void updateGP(base::DataVector knew, base::DataVector y);
	double acquisitionPI(base::DataVector knew, double kself, double bestsofar);
	double acquisitionEI(base::DataVector knew, double kself, double bestsofar);
  void CholeskyDecomposition();
  void solveCholeskySystem(base::DataVector& x);



    // double expkernel(base::DataVector x1, base::DataVector x2);
protected:
	base::DataMatrix kernelmatrix;
	base::DataMatrix kernelinv;
	base::DataVector transformedOutput;
	base::DataVector testknew;
	bool screwedvar;
	std::unique_ptr<optimization::FullSLE> sle;
	double maxofmax;
  base::DataMatrix gleft;
  // base::DataMatrix gright;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_ */
