/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DmModel.cpp
 *
 *  Created on: 15.04.2016
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/base/DmModule.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

using base::Grid;
using base::DataMatrix;
using base::DataVector;
using base::OperationEval;
using base::OperationMultipleEval;

DmModel::DmModel() : grid(nullptr), alpha(nullptr), lambda(0) {}

DmModel::~DmModel() {}

std::shared_ptr<Grid> DmModel::getGrid() { return this->grid; }

void sgpp::datadriven::DmModel::setGrid(std::shared_ptr<Grid> grid) { this->grid = grid; }

std::shared_ptr<DataVector> DmModel::getSurpluses() { return this->alpha; }

void sgpp::datadriven::DmModel::setSurplusses(std::shared_ptr<DataVector> alpha) {
  this->alpha = alpha;
}

double DmModel::getRegularizationParam() { return lambda; }

void DmModel::setRegularizationParam(double lambda) { this->lambda = lambda; }

double DmModel::evaluate(const DataVector& sample) {
  auto opEval(op_factory::createOperationEval(*grid));
  return opEval->eval(*alpha, sample);
}

std::unique_ptr<DataVector> DmModel::evaluate(DataMatrix& samples) {
  auto opMultEval(op_factory::createOperationMultipleEval(*grid, samples));
  auto result = std::make_unique<DataVector>(samples.getNrows());
  opMultEval->eval(*alpha, *result);
  return result;
}

double DmModel::score(Metric& scorer, DataMatrix& samples, const DataVector& trueValues) {
  auto predictedValues = evaluate(samples);
  return scorer(*predictedValues, trueValues);
}

} /* namespace datadriven */
} /* namespace sgpp */
