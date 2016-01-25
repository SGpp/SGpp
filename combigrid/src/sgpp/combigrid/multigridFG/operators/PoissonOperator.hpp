// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POISSONOPERATOR_HPP_
#define POISSONOPERATOR_HPP_

#include <combigrid.hpp>
#include <sgpp/combigrid/multigridFG/interface/OperatorFG.hpp>

namespace combigrid {

  /** abstract interface for the right hand side. */
  class CallBackRHS {
    public:
      /** empty CTor*/
      CallBackRHS() {
        ;
      }

      /** the callback function for the right hand side
       * @param coords the coordinates */
      virtual double eval(std::vector<double>&  coords) const {
        return 1.0;
      }
  };

  /** Class for constant right hand side*/
  class ConstRHS : public CallBackRHS {
    public:
      /** Ctor
       * @param v [IN] the constant value for the RHS*/
      ConstRHS(double v) : CallBackRHS() , constVal_(v) {
        ;
      }
      /** the callback function for the right hand side
       * @param coords [IN] the coordinates */
      virtual double eval(std::vector<double>&  coords) const {
        return constVal_;
      }
    private:
      /** const value for the right hand side*/
      double constVal_;
  };


  /** Poisson operator, solves the stationary Poisson equation with constant Dirichlet boundary condition */
  class PoissonOperator: public combigrid::OperatorFG {
    public:

      /** Ctor
       * @param fg [IN] the full grid on which this problem should be set up
       * @param sigma [IN] vector with the sigma values
       * @param callbackRHS [IN] the callback object for the right hand side*/
      PoissonOperator(const FullGridD* fg ,
                      const std::vector<double>& sigma ,
                      const CallBackRHS* callbackRHS );

      /** Dtor*/
      virtual ~PoissonOperator();

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


    private:

      /** lambda value used for regression*/
      std::vector<double> sigma_;

      /** dimension of the problem*/
      int dim_;

      /** */
      const CallBackRHS* callbackRHS_;

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

  };

}

#endif /* POISSONOPERATOR_HPP_ */