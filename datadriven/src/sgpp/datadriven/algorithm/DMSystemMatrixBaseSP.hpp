// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSPBASE_HPP
#define DMSYSTEMMATRIXSPBASE_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>
#include <sgpp/base/operation/hash/OperationMatrixSP.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

/**
 * Abstract class that defines the virtual class base::OperationMatrix for
 * classification and regression problems (single precision version)
 */
class DMSystemMatrixBaseSP : public base::OperationMatrixSP {
 protected:
  /// the dataset
  base::DataMatrixSP* dataset_;
  /// the lambda, the regularisation parameter
  float lambda_;
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
   * @param trainData matrix with training data
   * @param lambda the lambda, the regression parameter
   */
  DMSystemMatrixBaseSP(base::DataMatrixSP& trainData, float lambda);

  /**
   * Std-Destructor
   */
  virtual ~DMSystemMatrixBaseSP();

  virtual void mult(base::DataVectorSP& alpha,
                    base::DataVectorSP& result) = 0;

  /**
   * Generates the right hand side of the classification equation
   *
   * @param classes the class information of the training data
   * @param b reference to the vector that will contain the result of the matrix vector
   *   multiplication on the rhs
   */
  virtual void generateb(base::DataVectorSP& classes,
                         base::DataVectorSP& b) = 0;

  /**
   * forward declaration
   *
   * rebuilds the base::DataMatrix for Level and Index
   * this routine is needed for supporting adaptiva grids
   * with vectorized high performance kernels
   */
  virtual void rebuildLevelAndIndex();

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
  virtual void getTimers(double& timeMult, double& computeMult,
                         double& timeMultTrans, double& computeMultTrans);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* DMSYSTEMMATRIXSPBASE_HPP */
