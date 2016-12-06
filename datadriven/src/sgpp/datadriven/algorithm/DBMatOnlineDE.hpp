// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef DBMATONLINEDE_HPP_
#define DBMATONLINEDE_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>

#include <list>

/**
 * Class that captures a density function as an online object
 */
class DBMatOnlineDE : public DBMatOnline {
 public:
  /**
   * Constructor
   *
   * @param beta The initial weighting factor
   */

  explicit DBMatOnlineDE(double beta = 0.);
  /**
   * Destructor
   */
  virtual ~DBMatOnlineDE();

  /**
   * Reads an offline object
   *
   * @param o the offline object
   */
  virtual void readOffline(DBMatOffline* o);

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
  virtual void computeDensityFunction(
      sgpp::base::DataMatrix& m, bool save_b = false, bool do_cv = false,
      std::list<size_t>* deletedPoints = nullptr, size_t newPoints = 0);

  /**
   * Evaluates the density function at a certain point
   *
   * @param p the point at which the function is evaluated
   * @param force if set, it will even try to evaluate if the internal state
   * recommends otherwise
   */
  virtual double eval(sgpp::base::DataVector& p, bool force = false);

  /**
   * Return a pointer to alpha
   */
  virtual sgpp::base::DataVector* getAlpha();

  /**
   * Update alpha vector, i.e. delete entries specified by 'deletedPoints'
   * and/or
   * add 'newPoints' new entries
   *
   * @param deletedPoints indices of entries corresponding to deleted grid
   * points
   * @param newPoints number of added grid points
   */
  virtual void updateAlpha(std::list<size_t>* deletedPoints, size_t newPoints);

  /**
   * Returns if the surplus has already been computed
   */
  virtual bool isComputed();

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
  virtual void setCrossValidationParameters(int lambda_step,
                                            double lambda_start,
                                            double lambda_end,
                                            sgpp::base::DataMatrix* test,
                                            sgpp::base::DataMatrix* test_cc,
                                            bool logscale);

  /**
   * Returns the last best lamda
   */
  virtual double getBestLambda();

  /**
   * Sets the weighting factor
   *
   * @param beta the new weighting factor. If set to 0, no plasticity takes
   * place.
   */
  virtual void setBeta(double beta);

  /**
   * Returns the current weighting factor
   */
  virtual double getBeta();

  /**
   * Normalize the Density
   */
  virtual double normalize(size_t samples = 1000);

 protected:
  sgpp::base::DataVector* alpha_;
  bool functionComputed_;
  sgpp::base::DataVector* b_save;
  sgpp::base::DataVector* b_totalPoints;
  double beta_;
  size_t total_points;
  bool can_cv;
  int lambda_step_;
  double lambda_start_, lambda_end_;
  sgpp::base::DataMatrix *test_mat, *test_mat_res;
  double lambda;
  bool cv_logscale;
  double normFactor;
  size_t o_dim;

 private:
  double computeL2Error();
  double resDensity(sgpp::base::DataVector*& alpha);
};

#endif /* DBMATONLINEDE_HPP_ */

#endif /* USE_GSL */
