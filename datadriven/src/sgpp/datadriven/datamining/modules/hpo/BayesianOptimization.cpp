/*
 * BayesianOptimization.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include "BayesianOptimization.hpp"
#include "sgpp/optimization/sle/system/FullSLE.hpp"

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
	kernelinv.mult(knew, tmp);
	return kself - knew.dotProduct(tmp);
}

double BayesianOptimization::acquisitionEI(base::DataVector knew, double kself, double bestsofar){
    return (mean(knew) - bestsofar)/var(knew, kself);
}

void BayesianOptimization::updateGP(base::DataVector knew, base::DataVector y){
	size_t idx = kernelmatrix.appendRow();
	kernelmatrix.appendCol(knew);
	kernelmatrix.setRow(idx, knew);
    kernelinv.resize(idx+1,idx+1);
    optimization::FullSLE sle(kernelmatrix);
    optimization::sle_solver::Eigen solver{};
    base::DataMatrix identity(idx+1, idx+1, 0);
    for(size_t i=0; i<idx+1;i++){
        identity.set(i,i,1);
    }
    solver.solve(sle, identity, kernelinv);
	kernelinv.mult(y, transformedOutput);
}

} /* namespace datadriven */
} /* namespace sgpp */
