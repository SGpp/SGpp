// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>

#include <list>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
 * Class that stores, generates and manipulates a density function during online phase in on/off
 * learning.
 */
class DBMatOnlineDE : public DBMatOnline {
 public:
  /**
   * Constructor
   *
   * @param offline The offline object we base our evaluations on.
   * @param beta The initial weighting factor
   */
  explicit DBMatOnlineDE(DBMatOffline& offline, double beta = 0.);

  /**
   * Computes the density function for a certain data matrix
   *
   * @param m the matrix that contains the data points
   * @param save_b Indicates whether the old right hand side should be saved and
   *        combined with the new right hand side (aka streaming)
   * @param do_cv Indicates whether crossvalidation should take place
   * @param deletedPoints indicates the indices of removed grid points due to
   * coarsening
   * @param newPoints indicates the amount of added points due to refinement
   */
  void computeDensityFunction(DataMatrix& m, bool save_b = false, bool do_cv = false,
                              std::list<size_t>* deletedPoints = nullptr, size_t newPoints = 0);

  /**
   * Evaluates the density function at a certain point
   *
   * @param p the point at which the function is evaluated
   * @param force if set, it will even try to evaluate if the internal state recommends otherwise
   * @return the result of the evaluation
   */
  double eval(const DataVector& p, bool force = false);

  /**
   * Evaluates the density function on multiple points
   *
   * @param values the points at which the function is evaluated
   * @param results the result of the evaluation
   * @param force if set, it will even try to evaluate if the internal state recommends otherwise
   */
  void eval(DataMatrix& values, DataVector& results, bool force = false);

  /**
   * Return a reference to alpha
   *
   */
  DataVector& getAlpha();

  /**
   * Update alpha vector, i.e. delete entries specified by 'deletedPoints'
   * and/or
   * add 'newPoints' new entries
   *
   * @param deletedPoints indices of entries corresponding to deleted grid
   * points
   * @param newPoints number of added grid points
   */
  void updateAlpha(std::list<size_t>* deletedPoints, size_t newPoints);

  /**
   * Returns if the surplus has already been computed
   */
  bool isComputed();

  /**
   * Sets the crossvalidation parameters
   *
   * @param lambda_step how many different lambdas are tried out
   * @param lambda_start The smallest possible lambda
   * @param lambda_end The biggest possible lambda
   * @param test The test matrix. If NULL, the value old is not overwritten.
   * @param test_cc The results of the points in the test matrix.
   *        If NULL, the value old is not overwritten.
   * @param logscale Indicates whether the values between lambda_start and
   *        lambda_end are searched using logscale or not
   */
  void setCrossValidationParameters(int lambda_step, double lambda_start, double lambda_end,
                                    DataMatrix* test, DataMatrix* test_cc, bool logscale);

  /**
   * Returns the last best lamda
   */
  double getBestLambda();

  /**
   * Sets the weighting factor
   *
   * @param beta the new weighting factor. If set to 0, no plasticity takes
   * place.
   */
  void setBeta(double beta);

  /**
   * Returns the current weighting factor
   */
  double getBeta();

  /**
   * Normalize the Density
   */
  double normalize(size_t samples = 1000);

  /**
   * Normalize the Density using Quadrature
   */
  double normalizeQuadrature();

 protected:
  virtual void solveSLE(DataVector& b, bool do_cv) = 0;
  double computeL2Error();
  double resDensity(DataVector& alpha);

  DataVector alpha;
  bool functionComputed;
  DataVector bSave;
  DataVector bTotalPoints;
  double beta;
  size_t totalPoints;
  bool canCV;
  int lambdaStep;
  double lambdaStart, lambdaEnd;
  DataMatrix *testMat, *testMatRes;
  double lambda;
  bool cvLogscale;
  double normFactor;
  size_t oDim;
};

}  // namespace datadriven
}  // namespace sgpp
