// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CONJUGATEGRADIENTS_HPP
#define CONJUGATEGRADIENTS_HPP

#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace solver {

class ConjugateGradients : public SLESolver {
 public:
  /**
   * Std-Constructor
   */
  ConjugateGradients(size_t imax, double epsilon);

  /**
   * Std-Destructor
   */
  virtual ~ConjugateGradients();

  virtual void solve(sgpp::base::OperationMatrix& SystemMatrix, sgpp::base::DataVector& alpha,
                     sgpp::base::DataVector& b, bool reuse = false, bool verbose = false,
                     double max_threshold = -1.0);

  // Define functions for observer pattern in python

  /**
   * function that signals the start of the CG method (used in python)
   */
  virtual void starting();

  /**
   * function that signals the start of the calculation of the CG method (used in python)
   */
  virtual void calcStarting();

  /**
   * function that signals that one iteration step of the CG method has been completed (used in
   * python)
   */
  virtual void iterationComplete();

  /**
   * function that signals the finish of the cg method (used in python)
   */
  virtual void complete();
};

}  // namespace solver
}  // namespace sgpp

#endif /* CONJUGATEGRADIENTS_HPP */
