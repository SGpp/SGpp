/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataSourceState.hpp
 *
 *  Created on: 25.05.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <string>

namespace sgpp {
namespace datadriven {

struct DataSourceStateConfig {
  std::string filePath = "";
  size_t batchSize = 0;
  size_t numBatches = 0;
};

class DataSourceState {
 public:
  DataSourceState() : filePath(""), batchSize(0), numBatches(0), currentIteration(0) {}
  DataSourceState(const DataSourceStateConfig& config)
      : filePath(config.filePath),
        batchSize(config.batchSize),
        numBatches(config.numBatches),
        currentIteration(0) {}
  virtual ~DataSourceState(){};

  size_t getBatchSize() const { return batchSize; }

  size_t getCurrentIteration() const { return currentIteration; }

  void incrementCurrentIteration() { currentIteration++; }

  const std::string& getFilePath() const { return filePath; }

  size_t getNumBatches() const { return numBatches; }

 private:
  std::string filePath;
  size_t batchSize;
  size_t numBatches;
  size_t currentIteration;
};

} /* namespace datadriven */
} /* namespace sgpp */
