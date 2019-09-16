// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorCrossValidation.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

DataShufflingFunctorCrossValidation::DataShufflingFunctorCrossValidation(
    const CrossvalidationConfiguration crossValidationConfig,
    DataShufflingFunctor* shuffling) :
    shuffling{shuffling}, crossValidationConfig{crossValidationConfig}, currentFold{0} { }

DataShufflingFunctor* DataShufflingFunctorCrossValidation::clone() const {
  return new DataShufflingFunctorCrossValidation{*this};
}

void DataShufflingFunctorCrossValidation::setFold(size_t fold) {
  if (fold < crossValidationConfig.kfold_) {
    currentFold = fold;
  } else {
    throw new base::data_exception("Fold idx not availible!");
  }
}

size_t DataShufflingFunctorCrossValidation::operator()(size_t idx, size_t numSamples) {
  size_t foldSize = getCurrentFoldSize(numSamples);
  size_t foldStart = (numSamples / crossValidationConfig.kfold_) * currentFold;
  if (idx < foldSize) {
    // Map idx from {0, ..., foldSize - 1} to the {foldStart, ..., foldStart + foldSize - 1}
    return (*shuffling)(foldStart + idx, numSamples);
  }
  if (idx < foldSize + foldStart) {
    // Map idx from {foldSize, ..., foldSize + foldStart - 1} to {0, ..., foldStart - 1}
    return (*shuffling)(idx - foldSize, numSamples);
  } else {
    // Map idx from {foldSize + foldStart, ..., numSamples - 1} to {foldStart + foldSize, ...,
    // numSamples - 1}
    return (*shuffling)(idx, numSamples);
  }
}

size_t DataShufflingFunctorCrossValidation::getCurrentFoldSize(size_t numSamples) {
  size_t foldSize = numSamples / crossValidationConfig.kfold_;
  if (currentFold == crossValidationConfig.kfold_ - 1) {
    // The last fold possibly is bigger
    foldSize += numSamples % crossValidationConfig.kfold_;
  }
  return foldSize;
}

} /* namespace datadriven */
} /* namespace sgpp */








