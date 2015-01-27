// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef CONVERTPREWAVELETTOLINEAR_HPP
#define CONVERTPREWAVELETTOLINEAR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {



    /**
     * Class that implements the transformation of a hierarchical prewavelet sparse grid to a
     * hierarchical linear sparse grid. Therefore the ()operator is implemented in order to use
     * the sweep algorithm for the grid traversal. Let the coefficients from the hat basis be
     * \f$ h_{l,i}\f$ and from the prewavelet basis \f$ u_{l,i} \f$. To calculate the surplusses,
     * temp values are needed:
     * \f[
     * (l,i)\neq G_{n}^{1}:t_{l,i}=-\frac{6}{10}u_{l,i\pm1}+t_{l+1,2i}
     * \f]
     * All temp values for levels greater than the maximal level of the grid are set to 0. The actual
     * transformation is calculated as follows:
     * \f{eqnarray*}{
     * h_{l,i}&=&u_{l,i}+\frac{1}{10}u_{l,i\pm2}+t_{l+1,2i}-\frac{1}{2}t_{l,i\pm1}\qquad\mbox{if \ensuremath{(l,i)} is an inner point}\\h_{l,i}&=&\frac{9}{10}u_{l,i}+\frac{1}{10}u_{l,i\pm2}+t_{l+1,2i}-\frac{1}{2}t_{l,i\pm1}\qquad\mbox{if \ensuremath{(l,i)} is at border}
     * \f}
     *
     *The picture depicts all needed variables in oder to perform the transformation:
     * \image html prewavelets_dehierarch.png "This picture shows all involved gridpoints (red crosses) and temp values (green circles) to calculate the new hierarchical coefficients (red arrows) and new temp values (green arrows)."
     */
    class ConvertPrewaveletToLinear {
      protected:
        typedef GridStorage::grid_iterator grid_iterator;
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        /// the grid object
        GridStorage* storage;



      public:
        /**
         * Constructor, must be bind to a grid
         *
         * @param storage the grid storage object of the the grid, on which the hierarchisation should be executed
         */
        ConvertPrewaveletToLinear(GridStorage* storage);

        /**
         * Destructor
         */
        ~ConvertPrewaveletToLinear();
        /**
         * Converts a given prewavelet base to a normal linear base.
         */
        void operator()(DataVector& source, DataVector& result,
                        grid_iterator& index, size_t dim);
    };

    // namespace detail

  } // namespace SGPP
}

#endif /* CONVERTPREWAVELETTOLINEAR_HPP */