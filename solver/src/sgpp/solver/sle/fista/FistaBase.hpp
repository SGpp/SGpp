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
  virtual void solve(base::OperationMultipleEval& op, base::DataVector& weights,
                     const base::DataVector& b, size_t maxIt, double treshold) = 0;
};

}  // namespace solver
}  // namepsace sgpp

#endif  // FISTABASE_HPP
