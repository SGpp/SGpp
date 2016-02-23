// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditionalKDE.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <vector>

namespace SGPP {
namespace datadriven {

OperationDensityConditionalKDE::OperationDensityConditionalKDE(GaussianKDE& kde) : kde(&kde) {}

OperationDensityConditionalKDE::~OperationDensityConditionalKDE() {}

// -------------------------------------------------------------------

void OperationDensityConditionalKDE::doConditional(size_t mdim, float_t xbar,
                                                   GaussianKDE& conditionalizedKDE) {
  throw base::algorithm_exception(
      "OperationDensityConditionalKDE::doConditional is not implemented");
}

void OperationDensityConditionalKDE::doConditional(std::vector<size_t>& mdims,
                                                   base::DataVector& xbar,
                                                   datadriven::GaussianKDE& conditionalizedKDE) {
  // compute the dimensions to conditionalize
  size_t ndim = kde->getDim();
  std::vector<size_t> condDims(ndim - 1);

  size_t idim = 0;
  size_t jdim = 0;
  size_t kdim = 0;

  while (idim < ndim) {
    jdim = 0;

    while (jdim < mdims.size() && idim != mdims[jdim]) {
      jdim++;
    }

    // check if current dim has been found in dims
    if (jdim == mdims.size()) {
      condDims[kdim] = idim;
      kdim++;
    }

    idim++;
  }

  // update the conditionalization factors
  size_t nsamples = kde->getNsamples();
  base::DataVector pcond(nsamples);
  kde->getConditionalizationFactor(pcond);
  kde->updateConditionalizationFactors(xbar, condDims, pcond);

  // marginalize the kde to the desired dimensions
  op_factory::createOperationDensityMarginalizeKDE(*kde)->margToDimXs(mdims, conditionalizedKDE);
  // set the conditionalization coefficients
  conditionalizedKDE.setConditionalizationFactor(pcond);
}

void OperationDensityConditionalKDE::condToDimX(size_t mdim, base::DataVector& xbar,
                                                datadriven::GaussianKDE& conditionalizedKDE) {
  // compute the dimensions to conditionalize over
  size_t ndim = kde->getDim();
  std::vector<size_t> condDims(ndim - 1);
  size_t jdim = 0;

  for (size_t idim = 0; idim < ndim; idim++) {
    if (idim != mdim) {
      condDims[jdim] = idim;
      jdim++;
    }
  }

  // update the conditionalization factors
  size_t nsamples = kde->getNsamples();
  base::DataVector pcond(nsamples);
  kde->getConditionalizationFactor(pcond);
  kde->updateConditionalizationFactors(xbar, condDims, pcond);

  // marginalize the kde to the desired dimensions
  op_factory::createOperationDensityMarginalizeKDE(*kde)->margToDimX(mdim, conditionalizedKDE);
  // set the conditionalization coefficients
  conditionalizedKDE.setConditionalizationFactor(pcond);
}

void OperationDensityConditionalKDE::condToDimXs(std::vector<size_t>& mdims, base::DataVector& xbar,
                                                 datadriven::GaussianKDE& conditionalizedKDE) {
  throw base::algorithm_exception("OperationDensityConditionalKDE::condToDimXs is not implemented");
}

}  // namespace datadriven
}  // namespace SGPP
