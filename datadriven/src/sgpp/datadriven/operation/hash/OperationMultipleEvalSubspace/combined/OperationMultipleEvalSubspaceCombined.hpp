// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include "omp.h"
#include <immintrin.h>
#include <assert.h>

#include "OperationMultipleEvalSubspaceCombinedParameters.hpp"
#include "SubspaceNodeCombined.hpp"
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

class OperationMultipleEvalSubspaceCombined: public
  AbstractOperationMultipleEvalSubspace {
 private:

  sgpp::base::DataMatrix* paddedDataset;

  //size_t subspaceSize = -1;

  size_t maxGridPointsOnLevel;

  std::map<uint32_t, uint32_t> allLevelsIndexMap;

  size_t dim = -1;
  size_t maxLevel = 0;

  std::vector<SubspaceNodeCombined> allSubspaceNodes;
  uint32_t subspaceCount = -1;

  /// Pointer to the grid's gridstorage object
  //sgpp::base::GridStorage* storage = nullptr;
  uint32_t totalRegularGridPoints = -1;

#ifdef X86COMBINED_WRITE_STATS
  size_t refinementStep = 0;
  ofstream statsFile;
  string csvSep = "& ";
#endif

  void prepareSubspaceIterator();

  void listMultInner(size_t dim, const double* const datasetPtr,
                     sgpp::base::DataVector& alpha, size_t dataIndexBase,
                     size_t end_index_data, SubspaceNodeCombined& subspace,
                     double* levelArrayContinuous,
                     size_t validIndicesCount, size_t* validIndices, size_t* levelIndices,
                     //size_t *nextIterationToRecalcReferences, size_t nextIterationToRecalc,
                     double* evalIndexValuesAll, uint32_t* intermediatesAll);

  void uncachedMultTransposeInner(size_t dim, const double* const datasetPtr,
                                  size_t dataIndexBase,
                                  size_t end_index_data, SubspaceNodeCombined& subspace,
                                  double* levelArrayContinuous,
                                  size_t validIndicesCount, size_t* validIndices, size_t* levelIndices,
                                  //size_t *nextIterationToRecalcReferences,
                                  double* componentResults, double* evalIndexValuesAll,
                                  uint32_t* intermediatesAll);

  void setCoefficients(sgpp::base::DataVector& surplusVector);

  void unflatten(sgpp::base::DataVector& result);

  static uint32_t flattenIndex(size_t dim, std::vector<uint32_t>& maxIndices,
                               std::vector<uint32_t>& index);

  void setSurplus(std::vector<uint32_t>& level, std::vector<uint32_t>& maxIndices,
                  std::vector<uint32_t>& index,
                  double value);

  void getSurplus(std::vector<uint32_t>& level, std::vector<uint32_t>& maxIndices,
                  std::vector<uint32_t>& index,
                  double& value, bool& isVirtual);

  uint32_t flattenLevel(size_t dim, size_t maxLevel,
                        std::vector<uint32_t>& level);

 public:
#include "OperationMultipleEvalSubspaceCombined_calculateIndexCombined.hpp"

  OperationMultipleEvalSubspaceCombined(sgpp::base::Grid& grid,
                                        sgpp::base::DataMatrix& dataset);

  ~OperationMultipleEvalSubspaceCombined();

  void prepare() override;

  void multTransposeImpl(sgpp::base::DataVector& alpha,
                         sgpp::base::DataVector& result, const size_t start_index_data,
                         const size_t end_index_data) override;

  void multImpl(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                const size_t start_index_data,
                const size_t end_index_data) override;

  sgpp::base::DataMatrix* padDataset(sgpp::base::DataMatrix& dataset);

  size_t getAlignment() override;

  std::string getImplementationName() override;

  size_t getPaddedDatasetSize() override;
};

}
}
