// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/OperationPole.hpp>

namespace sgpp {
namespace combigrid {

class OperationPoleNodalisationBspline : public OperationPole {
 public:
  explicit OperationPoleNodalisationBspline(size_t degree, level_t level = 0,
      bool hasBoundary = true) : OperationPole(level, hasBoundary), degree(degree) {
  }

  ~OperationPoleNodalisationBspline() override {
  }

  void apply(base::DataVector& values, size_t start, size_t step, size_t count) override {
    switch (degree) {
      case 1: {
        // do nothing, as nodal coefficients equal values
        break;
      }
      case 3: {
        const double a = 1.0/6.0;
        const double b = 2.0/3.0;
        const double c = a;
        base::DataVector b2(count, b);
        base::DataVector d2(count);
        d2[0] = values[start];
        size_t j = start + step;

        for (size_t i = 1; i < count; i++) {
          const double w = a / b2[i-1];
          b2[i] -= w * c;
          d2[i] = values[j] - w * d2[i-1];
          j += step;
        }

        j -= step;
        values[j] = d2[count-1] / b2[count-1];

        for (size_t i = count-1; i-- > 0; ) {
          j -= step;
          values[j] = (d2[i] - c * values[j+step]) / b2[i];
        }

        break;
      }
    }
  }

 protected:
  size_t degree;
};

}  // namespace combigrid
}  // namespace sgpp
