// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef DBMATDECOMPMATRIXSOLVER_HPP_
#define DBMATDECOMPMATRIXSOLVER_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/SGSolver.hpp>

class DBMatDecompMatrixSolver : public sgpp::solver::SGSolver {
 public:
  DBMatDecompMatrixSolver();

  virtual ~DBMatDecompMatrixSolver();
};

#endif /* DBMATDECOMPMATRIXSOLVER_HPP_ */

#endif /* USE_GSL */
