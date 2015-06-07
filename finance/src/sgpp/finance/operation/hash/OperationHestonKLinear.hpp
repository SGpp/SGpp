// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHESTONKLINEAR_HPP
#define OPERATIONHESTONKLINEAR_HPP

#include <sgpp/pde/algorithm/UpDownFourOpDims.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    /**
     * Implements the Heston K-Operation (corresponds to operator K in Master's thesis), that is needed
     * the solve the multidimensional Heston
     * equation, on grids with fix Dirichlet-0-Boundaries.
     *
     */
    class OperationHestonKLinear : public SGPP::pde::UpDownFourOpDims {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         * @param coef vector that contains the constant coefficients of this operation
         */
        OperationHestonKLinear(SGPP::base::GridStorage* storage, float_t**** * coef);

        /**
         * Destructor
         */
        virtual ~OperationHestonKLinear();

      protected:

        /**
         * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the up-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the down-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * down-Gradient step in dimension <i>dim</i> applies the phi dphi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * up-Gradient step in dimension <i>dim</i> applies the phi dphi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * down-Gradient step in dimension <i>dim</i> applies the sqrt phi phi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * up-Gradient step in dimension <i>dim</i> applies the sqrt phi phi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * down-Gradient step in dimension <i>dim</i> applies the phi dphi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * up-Gradient step in dimension <i>dim</i> applies the phi dphi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * down-Gradient step in dimension <i>dim</i> applies the sqrt phi phi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * up-Gradient step in dimension <i>dim</i> applies the sqrt phi phi operation
         * in one dimension
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        /**
         * This method does nothing (this situation doesn't come up in Heston's PDEs). Needed only to make the class concrete.
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        void upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

    };

  }
}

#endif /* OPERATIONHESTONKLINEAR_HPP */