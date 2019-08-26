// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FISTABASE_HPP
#define FISTABASE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

namespace sgpp {
namespace solver {

class FistaBase {
 public:
  virtual ~FistaBase() {}
  virtual void solve(base::OperationMultipleEval& op, base::DataVector& weights,
                     const base::DataVector& classes, size_t maxIt, double treshold,
                     double L = 0.5) = 0;
  double getL() { return L; }

 protected:
  double L = 0.5;
};

}  // namespace solver
}  // namespace sgpp

#endif  // FISTABASE_HPP
