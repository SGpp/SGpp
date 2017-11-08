// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/SGSolver.hpp>

namespace sgpp {
namespace datadriven {

class DBMatDecompMatrixSolver : public sgpp::solver::SGSolver {
 public:
  DBMatDecompMatrixSolver();
};

}  // namespace datadriven
}  // namespace sgpp
