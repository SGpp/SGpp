// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NONSMOOTHFUNCTION_HPP
#define NONSMOOTHFUNCTION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

namespace sgpp {
namespace solver {

class RegularizationFunction {
 public:
  virtual double eval(base::DataVector weights) = 0;
  virtual base::DataVector prox(const base::DataVector& weights, double stepsize) = 0;
};

}  // namespace solver
}  // namespace sgpp

#endif  // NONSMOOTHFUNCTION_HPP
