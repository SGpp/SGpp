// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/tools/PartitioningTool.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class AbstractOperationMultipleEvalSubspace : public base::OperationMultipleEval {
 protected:
  base::GridStorage* storage;

 private:
  base::SGppStopwatch timer;
  double duration;

 public:
  AbstractOperationMultipleEvalSubspace(base::Grid& grid, base::DataMatrix& dataset)
      : base::OperationMultipleEval(grid, dataset), storage(&grid.getStorage()), duration(-1.0) {}

  ~AbstractOperationMultipleEvalSubspace() override {}

  virtual void multImpl(base::DataVector& alpha, base::DataVector& result,
                        const size_t start_index_data, const size_t end_index_data) = 0;

  virtual void multTransposeImpl(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                                 const size_t start_index_data, const size_t end_index_data) = 0;

  void multTranspose(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) override {
    if (!this->isPrepared) {
      this->prepare();
    }

    size_t originalAlphaSize = alpha.getSize();

    const size_t start_index_data = 0;
    const size_t end_index_data = this->getPaddedDatasetSize();

    // pad the alpha vector to the padded size of the dataset
    alpha.resizeZero(this->getPaddedDatasetSize());

    this->timer.start();
    result.setAll(0.0);

#pragma omp parallel
    {
      size_t start;
      size_t end;
      PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
                                                  this->getAlignment());
      this->multTransposeImpl(alpha, result, start, end);
    }

    alpha.resize(originalAlphaSize);
    this->duration = this->timer.stop();
  }

  void mult(sgpp::base::DataVector& source, sgpp::base::DataVector& result) override {
    if (!this->isPrepared) {
      this->prepare();
    }

    size_t originalResultSize = result.getSize();
    result.resizeZero(this->getPaddedDatasetSize());

    const size_t start_index_data = 0;
    const size_t end_index_data = this->getPaddedDatasetSize();

    this->timer.start();
    result.setAll(0.0);

#pragma omp parallel
    {
      size_t start;
      size_t end;
      PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
                                                  this->getAlignment());
      this->multImpl(source, result, start, end);
    }

    result.resize(originalResultSize);

    this->duration = this->timer.stop();
  }

  virtual size_t getPaddedDatasetSize() { return this->dataset.getNrows(); }

  virtual size_t getAlignment() = 0;

  virtual double getDuration() override { return this->duration; }

  static inline size_t getChunkGridPoints() { return 12; }
  static inline size_t getChunkDataPoints() {
    return 24;  // must be divisible by 24
  }
};
}
}
