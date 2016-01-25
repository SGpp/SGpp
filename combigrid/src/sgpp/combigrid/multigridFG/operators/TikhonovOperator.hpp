// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TIKHONOVOPERATOR_HPP_
#define TIKHONOVOPERATOR_HPP_

#include <combigrid.hpp>
#include <sgpp/combigrid/multigridFG/interface/OperatorFG.hpp>

namespace combigrid {

  /** The Tikhonov operator calss, which contains all the problem specific informations
   * and operations.*/
  class TikhonovOperator: public combigrid::OperatorFG {
    public:

      /** Ctor for a Tikhonov regularization problem (operator) on a given full grid
       * @param fg the full grid on which this problem should be set up
       * @param nrInputPoints nmber of input (e.g. Monte-Carlo) points
       * @param lambda lambda parameter
       * @param xCoords the coordinates of the input points for the regression (e.g. Monte-Carlo points). <br>
       *        the length of the vector is nrInputPoints*dimension , the small index is the dimension
       * @param yCoords the points values at the specified coordinates, the vector length is nrInputPoints */
      TikhonovOperator(const FullGridD* fg ,
                       int nrInputPoints ,
                       double lambda ,
                       const std::vector<double>* xCoords ,
                       const std::vector<double>* yCoords );

      /** Dtor*/
      virtual ~TikhonovOperator();

      /** method to create a new operator which acts on a given full grid
       * @param fg [IN] */
      virtual OperatorFG* factory(const FullGridD* fg) const;

      /** method to get the right hand side vector for one level
       * @param rhs [OUT] vector with the right hand side values
       * @param nrSpace [OUT] number of unknowns per node */
      virtual void getRHS(std::vector<double>& rhs , int& nrSpace) const;

      /** method to multiply one vector with the operator matrix.
       * @param inVect [IN] input vector
       * @param outVect [OUT] output vector (here will be the result A*inVect = outVect)*/
      virtual void multiplyVector(std::vector<double>& inVect , std::vector<double>& outVect) const;

      /** method to make Smoothing iterations (the number of iterations is specified by the user).
       * @param nrIt [IN] number of smoothing iterations to execute
       * @param u [IN/OUT] the unknown vector
       * @param rhs [IN] the right hand side for the iterations */
      virtual void doSmoothing(int nrIt ,
                               std::vector<double>& u, std::vector<double>& rhs) const;

      /** reset the lambda parameter */
      void setNewLambda(double lambda);

    private:

      /** lambda value used for regression*/
      double lambda_;

      /** dimension of the problem*/
      int dim_;

      // --------- CRS storage ----------
      /** number of rows an columns */
      int nrElem_;

      /** number of possible non-zero matrix element*/
      int nrMatrixElem_;

      /** CRS matrix storage: the matrix values*/
      std::vector<double> matrixVal_;

      /** CRS matrix storage: column index */
      std::vector<int> col_ind_;

      /** CRS matrix storage: row index */
      std::vector<int> row_ptr_;

      // ------------ the RHS (right hand side) ----
      std::vector<double> rhs_;

      // -------- the regression -------
      /** number of regression points */
      int nrRegPoints_;

      /** regression points */
      const std::vector<double>* xCoords_;
      /** values at the regression points */
      const std::vector<double>* yCoords_;
  };

}

#endif /* TIKHONOVOPERATOR_HPP_ */