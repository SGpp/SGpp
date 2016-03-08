// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef QUADRULE1D_HPP_
#define QUADRULE1D_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

class QuadRule1D {
 public:
  QuadRule1D();
  virtual ~QuadRule1D();

  void getLevelPointsAndWeights(size_t level, base::DataVector& coordinates,
                                base::DataVector& weights);

 protected:
  const size_t maxSupportedLevel = 20;
  std::vector<DataVector*> coordinates;
  std::vector<DataVector*> weights;
};

}  // namespace base
}  // namespace sgpp

#endif /* QUADRULE1D_HPP_ */
