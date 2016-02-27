// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CONJUGATEGRADIENTSMPI_HPP
#define CONJUGATEGRADIENTSMPI_HPP

#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace SGPP {
namespace parallel {

class ConjugateGradientsMPI : public SGPP::solver::SLESolver {
 private:
  /**
   * Routine called by the MPI slaves, here just the execution of
   * of sub part of the SystemMatrix's mult-Routine is needed.
   *
   * @param SystemMatrix reference to an SGPP::base::OperationMatrix Object that implements the
   * matrix vector multiplication
   * @param alpha the sparse grid's coefficients which have to be determined
   */
  virtual void waitForTask(SGPP::base::OperationMatrix& SystemMatrix,
                           SGPP::base::DataVector& alpha);

 public:
  explicit ConjugateGradientsMPI(size_t imax = 0, double epsilon = 0.0);

  virtual ~ConjugateGradientsMPI();

  virtual void solve(SGPP::base::OperationMatrix& SystemMatrix, SGPP::base::DataVector& alpha,
                     SGPP::base::DataVector& b, bool reuse = false, bool verbose = false,
                     double max_threshold = -1.0);
};
}  // namespace parallel
}  // namespace SGPP

#endif /* CONJUGATEGRADIENTSMPI_HPP */
