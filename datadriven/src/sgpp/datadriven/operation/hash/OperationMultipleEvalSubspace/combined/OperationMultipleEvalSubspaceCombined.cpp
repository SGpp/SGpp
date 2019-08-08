// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>

#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

OperationMultipleEvalSubspaceCombined::OperationMultipleEvalSubspaceCombined(Grid& grid,
                                                                             DataMatrix& dataset)
    : AbstractOperationMultipleEvalSubspace(grid, dataset) {
  this->paddedDataset = this->padDataset(dataset);
  this->storage = &grid.getStorage();
  // this->dataset = dataset;
  this->dim = dataset.getNcols();

  // this->subspaceSize = (2 * this->dim) + 8;
  this->maxGridPointsOnLevel = 0;

#ifdef X86COMBINED_WRITE_STATS
  string prefix("results/data/stats_");
  string fileName(X86COMBINED_WRITE_STATS);
  this->statsFile.open(prefix + fileName, ios::out);

  this->statsFile << "# name: " << X86COMBINED_WRITE_STATS_NAME << endl;
  this->statsFile << "refinementStep & ";
  this->statsFile << "nonVirtualGridPoints & ";
  this->statsFile << "totalRegularGridPoints & ";
  this->statsFile << "actualGridPoints & ";
  this->statsFile << "largestArraySubspace & ";
  this->statsFile << "largestListSubspace & ";
  this->statsFile << "numberOfListSubspaces & ";
  this->statsFile << "subspaceCount & ";
  this->statsFile << "avrPointsPerSubspace & ";
  this->statsFile << "memoryEstimate & ";
  this->statsFile << "memoryEfficiency";
  this->statsFile << endl;

#endif
}

OperationMultipleEvalSubspaceCombined::~OperationMultipleEvalSubspaceCombined() {
#ifdef X86COMBINED_WRITE_STATS
  this->statsFile.close();
#endif
}

void OperationMultipleEvalSubspaceCombined::prepare() {
  this->allLevelsIndexMap.clear();
  this->allSubspaceNodes.clear();

  this->prepareSubspaceIterator();
}

void OperationMultipleEvalSubspaceCombined::setCoefficients(DataVector& surplusVector) {
  std::vector<uint32_t> level(dim);
  std::vector<uint32_t> maxIndex(dim);
  std::vector<uint32_t> index(dim);

  base::level_t curLevel;
  base::index_t curIndex;

  for (size_t gridPoint = 0; gridPoint < this->storage->getSize(); gridPoint++) {
    sgpp::base::GridPoint& point = this->storage->getPoint(gridPoint);

    for (size_t d = 0; d < this->dim; d++) {
      point.get(d, curLevel, curIndex);
      level[d] = curLevel;
      index[d] = curIndex;
      maxIndex[d] = 1 << curLevel;
    }

    this->setSurplus(level, maxIndex, index, surplusVector.get(gridPoint));
  }
}

// writes a result vector in the order of the points in the grid storage
void OperationMultipleEvalSubspaceCombined::unflatten(DataVector& result) {
  std::vector<uint32_t> level(dim);
  std::vector<uint32_t> maxIndex(dim);
  std::vector<uint32_t> index(dim);

  base::level_t curLevel;
  base::index_t curIndex;

  for (size_t gridPoint = 0; gridPoint < this->storage->getSize(); gridPoint++) {
    sgpp::base::GridPoint& point = this->storage->getPoint(gridPoint);

    for (size_t d = 0; d < this->dim; d++) {
      point.get(d, curLevel, curIndex);
      level[d] = curLevel;
      index[d] = curIndex;
      maxIndex[d] = 1 << curLevel;
    }

    double surplus;
    bool isVirtual;
    this->getSurplus(level, maxIndex, index, surplus, isVirtual);

    result.set(gridPoint, surplus);
  }
}

void OperationMultipleEvalSubspaceCombined::setSurplus(std::vector<uint32_t>& level,
                                                       std::vector<uint32_t>& maxIndices,
                                                       std::vector<uint32_t>& index, double value) {
  uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
  uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
  uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
  SubspaceNodeCombined& subspace = this->allSubspaceNodes[subspaceIndex];
  subspace.setSurplus(indexFlat, value);
}

void OperationMultipleEvalSubspaceCombined::getSurplus(std::vector<uint32_t>& level,
                                                       std::vector<uint32_t>& maxIndices,
                                                       std::vector<uint32_t>& index, double& value,
                                                       bool& isVirtual) {
  uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
  uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
  uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
  SubspaceNodeCombined& subspace = this->allSubspaceNodes[subspaceIndex];
  value = subspace.getSurplus(indexFlat);

  if (std::isnan(value)) {
    isVirtual = true;
  } else {
    isVirtual = false;
  }
}

uint32_t OperationMultipleEvalSubspaceCombined::flattenIndex(size_t dim,
                                                             std::vector<uint32_t>& maxIndices,
                                                             std::vector<uint32_t>& index) {
  uint32_t indexFlat = index[0];
  indexFlat >>= 1;

  for (size_t i = 1; i < dim; i++) {
    uint32_t actualDirectionGridPoints = maxIndices[i];
    actualDirectionGridPoints >>= 1;
    indexFlat *= actualDirectionGridPoints;
    uint32_t actualIndex = index[i];
    actualIndex >>= 1;  // divide index by 2, skip even indices
    indexFlat += actualIndex;
  }

  return indexFlat;
}

uint32_t OperationMultipleEvalSubspaceCombined::flattenLevel(size_t dim, size_t maxLevel,
                                                             std::vector<uint32_t>& level) {
  uint32_t levelFlat = 0;
  levelFlat += level[dim - 1];

  // loop terminates at -1
  for (int i = static_cast<int>(dim - 2); i >= 0; i--) {
    levelFlat *= static_cast<uint32_t>(maxLevel);
    levelFlat += level[i];
  }

  return levelFlat;
}

DataMatrix* OperationMultipleEvalSubspaceCombined::padDataset(sgpp::base::DataMatrix& dataset) {
  size_t chunkSize = X86COMBINED_PARALLEL_DATA_POINTS;

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = dataset.getNrows() % chunkSize;
  size_t loopCount = chunkSize - remainder;

  if (loopCount == chunkSize) {
    return &dataset;
  }

  sgpp::base::DataVector lastRow(dataset.getNcols());
  size_t oldSize = dataset.getNrows();
  dataset.getRow(dataset.getNrows() - 1, lastRow);

  DataMatrix* paddedDataset = new DataMatrix(dataset);
  // pad to make: dataset % X86COMBINED_PARALLEL_DATA_POINTS == 0
  paddedDataset->resize(dataset.getNrows() + loopCount);

  for (size_t i = 0; i < loopCount; i++) {
    paddedDataset->setRow(oldSize + i, lastRow);
  }

  // additional padding for subspace skipping
  // if validIndices contain X86COMBINED_PARALLEL_DATA_POINTS - 1 it is possible for a vector
  // iteration to contain
  // indices larger than size(dataset) (even though the dataset is divided by
  // X86COMBINED_PARALLEL_DATA_POINTS)
  // add X86COMBINED_VEC_PADDING dummy data points to avoid that problem
  // add X86COMBINED_VEC_PADDING * 2 to also enable the calculateIndexCombined2() method
  // this works due to special semantics of "reserveAdditionalRows()", this function adds additional
  // unused (and uncounted) rows
  paddedDataset->reserveAdditionalRows(X86COMBINED_VEC_PADDING * 2);

  for (size_t i = paddedDataset->getNrows();
       i < paddedDataset->getNrows() + paddedDataset->getAdditionallyReservedRows(); i++) {
    for (size_t j = 0; j < paddedDataset->getNcols(); j++) {
      paddedDataset->set(i, j, 0.0);
    }
  }

  return paddedDataset;
}

size_t OperationMultipleEvalSubspaceCombined::getPaddedDatasetSize() {
  return this->paddedDataset->getNrows();
}

size_t OperationMultipleEvalSubspaceCombined::getAlignment() {
  return X86COMBINED_PARALLEL_DATA_POINTS;
}

std::string OperationMultipleEvalSubspaceCombined::getImplementationName() { return "COMBINED"; }

}  // namespace datadriven
}  // namespace sgpp
