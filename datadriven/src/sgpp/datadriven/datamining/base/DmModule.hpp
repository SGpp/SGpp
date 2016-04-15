/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Module.hpp
 *
 *  Created on: Apr 12, 2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <memory>
#include <sgpp/datadriven/datamining/base/DmStateStorage.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

namespace sgpp {
namespace datadriven {

class DmModule {
 public:
  DmModule() : nextModule(nullptr), stateStorage(nullptr) {}
  DmModule(std::shared_ptr<DmModule> newNextModule, std::shared_ptr<DmStateStorage> sharedState)
      : nextModule(newNextModule), stateStorage(sharedState) {}
  virtual ~DmModule() {}
  virtual void run() = 0;
  void setNextModule(std::shared_ptr<DmModule> newNextModule) { nextModule = newNextModule; }
  std::shared_ptr<DmModule> getNextModule() { return nextModule; }

  std::shared_ptr<DmStateStorage> getStateStorage() { return stateStorage; }

  void setStateStorage(std::shared_ptr<DmStateStorage>& stateStorage) {
    this->stateStorage = stateStorage;
  }

 protected:
  std::shared_ptr<DmModule> nextModule;
  std::shared_ptr<DmStateStorage> stateStorage;
};
} /* namespace datadriven */
} /* namespace sgpp */
