/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DmConfig.cpp
 *
 *  Created on: 15.04.2016
 *      Author: Michael Lettrich
 */

#include "DmStateStorage.hpp"

namespace sgpp {
namespace datadriven {

DmStateStorage::~DmStateStorage() {
  // TODO(lettrich): Auto-generated destructor stub
}

DmStateStorage::DmStateStorage() {
  // TODO(lettrich): Auto-generated constructor stub
}

std::shared_ptr<Dataset> DmStateStorage::getDataset() { return dataset; }

void DmStateStorage::setDataset(std::shared_ptr<Dataset> dataset) { this->dataset = dataset; }

std::shared_ptr<DmModel> DmStateStorage::getModel() { return model; }

} /* namespace datadriven */
} /* namespace sgpp */
