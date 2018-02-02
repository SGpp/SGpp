/*
 * BayesianOptimization.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#include "BayesianOptimization.hpp"

namespace sgpp {
namespace datadriven {

BayesianOptimization::BayesianOptimization() {
	// TODO Auto-generated constructor stub

}

double BayesianOptimization::mean(base::DataVector knew){
	return knew.dotProduct(transformedOutput);
}

double BayesianOptimization::var(base::DataVector knew, double kself){
	base::DataVector tmp(knew.size());
	kernelmatrix.mult(knew, tmp);
	return kself - knew.dotProduct(tmp);
}

void BayesianOptimization::updateGP(base::DataVector knew, base::DataVector y){
	size_t idx = kernelmatrix.appendRow();
	kernelmatrix.appendCol(knew);
	kernelmatrix.setRow(idx, knew);
	//EDIT: inverse of kernelmatrix, two matrices needed
	kernelmatrix.mult(y, transformedOutput);
}

} /* namespace datadriven */
} /* namespace sgpp */
