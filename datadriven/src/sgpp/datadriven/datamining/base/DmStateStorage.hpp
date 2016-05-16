/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DmConfig.hpp
 *
 *  Created on: 15.04.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <memory>
#include <sgpp/datadriven/datamining/base/DmModel.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

namespace sgpp {
namespace datadriven {

class DmStateStorage {
 public:
  virtual ~DmStateStorage();
  DmStateStorage();
  std::shared_ptr<Dataset> getDataset();
  void setDataset(std::shared_ptr<Dataset> dataset);
  DmModel& getModel();

 private:
  std::shared_ptr<Dataset> dataset;
  DmModel model;
};

} /* namespace datadriven */
} /* namespace sgpp */
