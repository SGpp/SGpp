// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace combigrid {

class GeneralOperation {
 public:
  GeneralOperation() {
  }

  virtual ~GeneralOperation() {
  }

  virtual void apply(base::DataVector& values) = 0;

  virtual void apply(base::DataVector& values, size_t start, size_t step, size_t stop) {
    base::DataVector valuesCopy((stop - start - 1) / step + 1);
    size_t j = start;

    for (size_t i = 0; i < valuesCopy.size(); i++) {
      valuesCopy[i] = values[j];
      j += step;
    }

    apply(valuesCopy);

    j = start;

    for (size_t i = 0; i < valuesCopy.size(); i++) {
      values[j] = valuesCopy[i];
      j += step;
    }
  }
};

}  // namespace combigrid
}  // namespace sgpp
