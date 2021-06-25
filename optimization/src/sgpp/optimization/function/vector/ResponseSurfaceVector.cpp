// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/function/vector/ResponseSurfaceVector.hpp>

#include <string>

namespace sgpp {
namespace optimization {

sgpp::base::DataVector ResponseSurfaceVector::eval(sgpp::base::DataVector v) {
  sgpp::base::DataVector evaluations(numRes);
  interpolants->eval(v, evaluations);
  return evaluations;
}

// sgpp::base::DataVector ResponseSurfaceVector::evalJacobian(sgpp::base::DataVector v,
//                                                            sgpp::base::DataMatrix& jacobian) {
//   sgpp::base::DataVector evaluations(numRes);
//   interpolantGradients->eval(v, evaluations, jacobian);
//   return evaluations;
// }

sgpp::base::DataVector ResponseSurfaceVector::averageL2Error(
    std::shared_ptr<sgpp::base::VectorFunction> objectiveFunc,
    sgpp::base::DataVector& componentwiseErrors, size_t numMCPoints) {
  componentwiseErrors.resizeZero(numRes);
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] = sgpp::base::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
    }
    sgpp::base::DataVector evalInterpolant = this->eval(randomVector);
    sgpp::base::DataVector evalObjectiveFunc(numRes);
    objectiveFunc->eval(randomVector, evalObjectiveFunc);
    // componentwise l2 errors (f-tilde f)^2
    evalInterpolant.sub(evalObjectiveFunc);
    evalInterpolant.componentwise_mult(evalInterpolant);
    componentwiseErrors.add(evalInterpolant);
  }
  double averageL2Err = 0;
  for (size_t r = 0; r < numRes; r++) {
    componentwiseErrors[r] = sqrt(componentwiseErrors[r] / static_cast<double>(numMCPoints));
    averageL2Err += fabs(componentwiseErrors[r]);
  }
  averageL2Err /= static_cast<double>(numRes);
  // find minimum and maximum error over all outputs
  sgpp::base::DataVector compwiseErrCopy(componentwiseErrors);
  compwiseErrCopy.abs();
  sgpp::base::DataVector returnVec(3);
  returnVec[0] = averageL2Err;
  returnVec[1] = compwiseErrCopy.min();
  returnVec[2] = compwiseErrCopy.max();
  return returnVec;
}

sgpp::base::DataVector ResponseSurfaceVector::averageNRMSE(
    std::shared_ptr<sgpp::base::VectorFunction> objectiveFunc,
    sgpp::base::DataMatrix& componentwiseErrorData, size_t numMCPoints) {
  sgpp::base::DataVector componentwiseL2Errors(numRes);
  sgpp::base::DataVector componentwiseNRMSE(numRes);
  sgpp::base::DataVector minEvaluations(numRes, std::numeric_limits<double>::max());
  sgpp::base::DataVector maxEvaluations(numRes, std::numeric_limits<double>::lowest());
  sgpp::base::DataVector randomVector(numDim);
  for (size_t i = 0; i < numMCPoints; i++) {
    for (size_t d = 0; d < numDim; d++) {
      randomVector[d] = sgpp::base::RandomNumberGenerator::getInstance().getUniformRN(lb[d], ub[d]);
    }
    sgpp::base::DataVector evalInterpolant = this->eval(randomVector);
    sgpp::base::DataVector evalObjectiveFunc(numRes);
    objectiveFunc->eval(randomVector, evalObjectiveFunc);

    // todo(rehmemk) This seems like a very bad implementation. Improve it!
    for (size_t t = 0; t < numRes; t++) {
      if (evalObjectiveFunc[t] > maxEvaluations[t]) {
        maxEvaluations[t] = evalObjectiveFunc[t];
      }
      if (evalObjectiveFunc[t] < minEvaluations[t]) {
        minEvaluations[t] = evalObjectiveFunc[t];
      }
    }

    // componentwise l2 errors (f-tilde f)^2
    evalInterpolant.sub(evalObjectiveFunc);
    evalInterpolant.componentwise_mult(evalInterpolant);
    componentwiseL2Errors.add(evalInterpolant);
  }
  double averageL2Err = 0;
  double averageNRMSE = 0;
  for (size_t t = 0; t < numRes; t++) {
    componentwiseL2Errors[t] = sqrt(componentwiseL2Errors[t] / static_cast<double>(numMCPoints));
    if ((maxEvaluations[t] - minEvaluations[t]) > 1e-14) {
      componentwiseNRMSE[t] = componentwiseL2Errors[t] / (maxEvaluations[t] - minEvaluations[t]);
    } else {
      componentwiseNRMSE[t] = componentwiseL2Errors[t];
    }
    averageL2Err += componentwiseL2Errors[t];
    averageNRMSE += componentwiseNRMSE[t];
  }
  averageL2Err /= static_cast<double>(numRes);
  averageNRMSE /= static_cast<double>(numRes);

  componentwiseErrorData.resize(4, numRes);
  componentwiseErrorData.setRow(0, componentwiseNRMSE);
  componentwiseErrorData.setRow(1, componentwiseL2Errors);
  componentwiseErrorData.setRow(2, minEvaluations);
  componentwiseErrorData.setRow(3, maxEvaluations);

  // find minimum and maximum error over all outputs
  sgpp::base::DataVector compwiseNRMSECopy(componentwiseNRMSE);
  compwiseNRMSECopy.abs();
  sgpp::base::DataVector returnVec(4);
  returnVec[0] = averageNRMSE;
  returnVec[1] = averageL2Err;
  returnVec[2] = compwiseNRMSECopy.min();
  returnVec[3] = compwiseNRMSECopy.max();
  return returnVec;
}

void ResponseSurfaceVector::transformPoint(sgpp::base::DataVector& v,
                                           sgpp::base::DataVector lBounds,
                                           sgpp::base::DataVector uBounds,
                                           sgpp::base::DataVector newlBounds,
                                           sgpp::base::DataVector newuBounds) {
  v.sub(lBounds);
  uBounds.sub(lBounds);
  v.componentwise_div(uBounds);
  newuBounds.sub(newlBounds);
  v.componentwise_mult(newuBounds);
  v.add(newlBounds);
}

double ResponseSurfaceVector::domainVolume() {
  double vol = 1;
  for (size_t d = 0; d < numDim; d++) {
    vol *= ub[d] - lb[d];
  }
  return vol;
}

}  // namespace optimization
}  // namespace sgpp
