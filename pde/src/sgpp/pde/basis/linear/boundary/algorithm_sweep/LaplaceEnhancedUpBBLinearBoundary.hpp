// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LAPLACEENHANCEDUPBBLINEARBOUNDARY_HPP
#define LAPLACEENHANCEDUPBBLINEARBOUNDARY_HPP

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedUpBBLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implementation of sweep operator () for
     * enhanced Laplace operator, up operation.
     *
     * This sweep operator calculates all ups (L2 scalar products and
     * gradient) for a given dimension.
     */
    class LaplaceEnhancedUpBBLinearBoundary : public LaplaceEnhancedUpBBLinear {
      private:
        /**
         * calculates the L2 up operation on level 0 basis functions
         *
         * @param fl source coefficient of the left boundary basis function
         * @param fr source coefficient of the right boundary basis function
         * @param seq_left unique order number of the left boundary basis function
         * @param seq_right unique order number of the right boundary basis function
         * @param dim current dimension
         * @param algo_dim current algorithmic dimension
         * @param q stretching of basis function in the current algorithmic dimension
         */
        void calcL2Boundary(float_t fl, float_t fr, size_t seq_left, size_t seq_right, size_t dim, size_t algo_dim, float_t q);

        /**
         * calculates the gradient up operation on level 0 basis functions
         *
         * @param fl source coefficient of the left boundary basis function
         * @param fr source coefficient of the right boundary basis function
         * @param seq_left unique order number of the left boundary basis function
         * @param seq_right unique order number of the right boundary basis function
         * @param dim current dimension
         * @param algo_dim current algorithmic dimension
         * @param q_reci reciprocal of stretching of basis function in the current algorithmic dimension
         */
        void calcGradBoundary(float_t fl, float_t fr, size_t seq_left, size_t seq_right, size_t dim, size_t algo_dim, float_t q_reci);

      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        LaplaceEnhancedUpBBLinearBoundary(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~LaplaceEnhancedUpBBLinearBoundary();

        /**
         * This operations performs the calculation of up in the direction of dimension <i>dim</i>
         * on a grid with fix Dirichlet 0 boundary conditions
         *
         * @param source SGPP::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result SGPP::base::DataVector that contains the result of the up operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(SGPP::base::DataMatrix& source, SGPP::base::DataMatrix& result, grid_iterator& index, size_t dim);
    };

    // namespace detail
  }
  // namespace SGPP
}

#endif /* LAPLACEENHANCEDUPBBLINEARBOUNDARY_HPP */
