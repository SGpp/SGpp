// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>

namespace sgpp {
namespace datadriven {

DBMatDecompMatrixSolver::DBMatDecompMatrixSolver() : SGSolver(0, 0) {}

DBMatDecompMatrixSolver::~DBMatDecompMatrixSolver() {}

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */
