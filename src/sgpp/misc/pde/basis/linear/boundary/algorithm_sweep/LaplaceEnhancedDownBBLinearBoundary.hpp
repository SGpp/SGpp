/* *****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LAPLACEENHANCEDDOWNBBLINEARBOUNDARY_HPP
#define LAPLACEENHANCEDDOWNBBLINEARBOUNDARY_HPP

#include "misc/pde/basis/linear/noboundary/algorithm_sweep/LaplaceEnhancedDownBBLinear.hpp"

namespace sg {
  namespace pde {

    /**
     * Implementation of sweep operator () for
     * enhanced Laplace operator, down operation.
     *
     * This sweep operator calculates all downs (L2 scalar products and
     * gradient) for a given dimension.
     */
    class LaplaceEnhancedDownBBLinearBoundary : public LaplaceEnhancedDownBBLinear {
      private:
        /**
         * calculates the L2 down operation on level 0 basis functions
         *
         * @param fl source coefficient of the left boundary basis function
         * @param fr source coefficient of the right boundary basis function
         * @param seq_left unique order number of the left boundary basis function
         * @param seq_right unique order number of the right boundary basis function
         * @param dim current dimension
         * @param algo_dim current algorithmic dimension
         * @param q stretching of basis function in the current algorithmic dimension
         */
        void calcL2Boundary(double fl, double fr, size_t seq_left, size_t seq_right, size_t dim, size_t algo_dim, double q);

        /**
         * calculates the gradient down operation on level 0 basis functions
         *
         * @param fl source coefficient of the left boundary basis function
         * @param fr source coefficient of the right boundary basis function
         * @param seq_left unique order number of the left boundary basis function
         * @param seq_right unique order number of the right boundary basis function
         * @param dim current dimension
         * @param algo_dim current algorithmic dimension
         * @param q_reci reciprocal of stretching of basis function in the current algorithmic dimension
         */
        void calcGradBoundary(double fl, double fr, size_t seq_left, size_t seq_right, size_t dim, size_t algo_dim, double q_reci);

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        LaplaceEnhancedDownBBLinearBoundary(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~LaplaceEnhancedDownBBLinearBoundary();

        /**
         * This operations performs the calculation of down in the direction of dimension <i>dim</i>
         * on a grid with Dirichlet 0 boundary conditions.
         *
         * @param source sg::base::DataMatrix that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result sg::base::DataMatrix that contains the result of the down operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction', here all downs are calculated
         */
        virtual void operator()(sg::base::DataMatrix& source, sg::base::DataMatrix& result, grid_iterator& index, size_t dim);
    };

    // namespace detail
  }
  // namespace sg
}

#endif /* LAPLACEENHANCEDDOWNBBLINEARBOUNDARY_HPP */
