/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONREGULARIZATIONDIAGONAL_HPP
#define OPERATIONREGULARIZATIONDIAGONAL_HPP

#include "base/operation/OperationMatrix.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace datadriven {

    /**
     * Implementation of the application of a diagonal matrix to a
     * DataVector for regularization.
     * This class implements several scaling possibilities.
     */
    class OperationRegularizationDiagonal: public base::OperationMatrix {
      protected:
        // to remember mode
        int mode;
        // to remember parameter k
        double k;
        // to remember state of grid in terms of number of grid points
        size_t size;
        // ro remember grid's storage
        base::GridStorage* storage;
        // to store diagonal
        base::DataVector diagonal;

        // Diagonal entries have been provided
        //static const int DIAGMATRIX = 0;

        /**
         * Check for mode and run corresponding init
         */
        virtual void init();

        /**
         * Initialize Hkmix.
         * Diagonal entries are
         * @f$\prod_{k=1}^d \langle \phi_k(x_k),\phi_k(x_k) \rangle_{H^k}@f$.
         * @param k Parameter k
         */
        virtual void initHkmix(double k) = 0;

        /**
         * Initialize H0HkLaplace.
         * Diagonal entries are
         * @f$\sum_{k=1}^d \langle \phi_k(x_k),\phi_k(x_k) \rangle_{H^k}
         * \prod_{l\neq k} \langle \phi_k(x_k),\phi_k(x_k) \rangle_{H^0}@f$.
         * @param k Parameter k
         */
        virtual void initH0HkLaplace(double k) = 0;

        /**
         * Initialize ISOTROPIC_PENALTY, ignores constructor parameter k.
         * Diagonal entries are
         * @f$\frac{1}{\max\{l_1,\dots,l_d\}-\min\{l_1,\dots,l_d\}+1}d@f$.
         */
        virtual void initIsotropicPenalty();

        /**
         * Initialize ANISOTROPIC_PENALTY, ignores constructor parameter k.
         * Diagonal entries are
         * @f$\frac{1}{2}\log(1+(\frac{\max\{l_1,\dots,l_d\}}{\max\{\min\{l_1,\dots,l_d\},1\}})d)@f$.
         */
        virtual void initAnisotropicPenalty();


      public:
        /// Diagonal scaling for @f$H^k_\textbf{mix}@f$-norm (product of
        /// @f$H^k@f$ 1D norms)
        static const int HKMIX = 1;
        /// Diagonal scaling for pseudo-laplace, replacing @f$L_2@f$ with
        /// @f$H^0@f$ in 1d
        static const int H0HKLAPLACE = 2;
        // Penalizes strongly isotropic subspaces more
        static const int ISOTROPIC_PENALTY = 3;
        // Penalizes strongly anisotropic subspaces more
        static const int ANISOTROPIC_PENALTY = 4;

        /**
         * Constructor of OperationRegularizationDiagonal.
         * Constructor should most likely call init() in subclasses.
         * @param storage Pointer to grid's storage object
         * @param mode Mode, specifying which regularization to use. Example: OperationRegularizationDiagonal::HKMIX.
         * @param k Parameter for @f$H^k@f$
         */
        OperationRegularizationDiagonal(base::GridStorage* storage, int mode, double k);

        /**
         * Destructor
         */
        virtual ~OperationRegularizationDiagonal() {};

        /**
         * Multiplication with diagonal matrix.
         * @param alpha Source vector
         * @param result Result vector (entries are overwritten)
         */
        virtual void mult(base::DataVector& alpha, base::DataVector& result);
    };

  }
}
#endif /* OPERATIONREGULARIZATIONDIAGONAL_HPP */
