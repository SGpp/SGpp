// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorSequential.hpp>

namespace sgpp {
namespace datadriven {

DataShufflingFunctor* DataShufflingFunctorSequential::clone() const {
  return new DataShufflingFunctorSequential{*this};
}

size_t DataShufflingFunctorSequential::operator()(size_t idx, size_t numSamples) { return idx; }

} /* namespace datadriven */
} /* namespace sgpp */




