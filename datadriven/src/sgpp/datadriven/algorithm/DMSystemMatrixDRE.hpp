// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Abstract class that defines the virtual class base::OperationMatrix for the density ratio
 * estimation problem
 */
class DMSystemMatrixDRE : public base::OperationMatrix {
 protected:
  /// the datasets
  base::DataMatrix datasetP_;
  base::DataMatrix datasetQ_;
  /// the lambda, the regularisation parameter
  double lambda_;
  /// time needed for Mult
  double completeTimeMult_;
  /// time needed only for the computation of mult, interesting on accelerator boards
  double computeTimeMult_;
  /// time needed for Mult transposed
  double completeTimeMultTrans_;
  /// time needed only for the computation of mult transposed, interesting on accelerator boards
  double computeTimeMultTrans_;
  /// Stopwatch needed to determine the durations of mult and mult transposed
  base::SGppStopwatch* myTimer_;

 public:
  /**
   * Std-Constructor
   *
   * @param trainDataP matrix with training data for first density
   * @param trainDataQ matrix with training data for second density
   * @param lambda the lambda, the regression parameter
   */
  DMSystemMatrixDRE(base::DataMatrix& trainDataP, base::DataMatrix& trainDataQ, double lambda);

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixDRE();

  virtual void mult(base::DataVector& alpha, base::DataVector& result) = 0;

  /**
   * Generates the right hand side of the equation
   *
   * @param b reference to the vector that will contain the result of the matrix vector
   * multiplication on the rhs
   */
  virtual void generateb(base::DataVector& b) = 0;

  /**
   * forward declaration
   *
   * rebuilds the base::DataMatrix for Level and Index
   * this routine is needed for supporting adaptive grids
   * with vectorized high performance kernels
   */
  virtual void prepareGrid();

  /**
   * resets all timers to 0
   */
  virtual void resetTimers();

  /**
   * gets the timer's values by saving them into call by reference values
   *
   * @param timeMult variable to store overall time needed for Mult
   * @param computeMult variable to store compute time needed for Mult
   * @param timeMultTrans variable to store everall time needed for Mult Transposed
   * @param computeMultTrans variable to store compute time needed for Mult Transposed
   */
  virtual void getTimers(double& timeMult, double& computeMult, double& timeMultTrans,
                         double& computeMultTrans);
};

}  // namespace datadriven
}  // namespace sgpp
