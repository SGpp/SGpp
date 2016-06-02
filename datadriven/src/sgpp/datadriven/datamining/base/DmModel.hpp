/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DmModel.hpp
 *
 *  Created on: 15.04.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <memory>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

using base::Grid;
using base::DataMatrix;
using base::DataVector;

class DmModel {
 public:
  DmModel();
  DmModel(DmModel&& model);
  virtual ~DmModel();

  std::shared_ptr<Grid> getGrid();
  void setGrid(std::shared_ptr<Grid> grid);
  std::shared_ptr<DataVector> getSurpluses();
  void setSurplusses(std::shared_ptr<DataVector> alpha);
  double getRegularizationParam();
  void setRegularizationParam(double lambda);
  double evaluate(const DataVector& sample);
  std::unique_ptr<DataVector> evaluate(DataMatrix& samples);
  double score(Metric& scorer, DataMatrix& samples, const DataVector& values);

 private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<DataVector> alpha;
  double lambda;
};

} /* namespace datadriven */
} /* namespace sgpp */
