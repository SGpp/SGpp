/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de)
// @author xander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LAPLACEUPGRADIENTPREWAVELET_HPP
#define LAPLACEUPGRADIENTPREWAVELET_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {



    /**
     * Implements the upGradient Method needed for the Laplace operator on prewavelet grids. The calculation
     * is done iterative and utilizes the following temp variables:
     * \f[
     * t_{k,j}=-\frac{6}{10}u_{k,j\pm1}+t_{k+1,2j}\qquad(k,j)\notin G_{n}^{1}
     * \f]
     * The correct values are then calculated as follows:
     * \f{eqnarray*}{
     * r_{k,j}&=&\frac{1}{h_{k}}\left(2t_{k+1,2j}-t_{k+1,2(j\pm1)}\right)\\&&-\frac{6}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm2)}+2t_{k+1,2(j\pm1)}-2t_{k+1,2j}\right)\\&&+\frac{1}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm1)}+2t_{k+1,2(j\pm2)}-t_{k+1,2(j\pm2)}\right)
     * \f}
     * In case of borders:
     * \f{eqnarray*}{
     * r_{k,j}&=&\frac{9}{10}\frac{1}{h_{k}}\left(2t_{k+1,2j}-t_{k+1,2(j\pm1)}\right)\\&&-\frac{6}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm2)}+2t_{k+1,2(j\pm1)}-t_{k+1,2j}\right)\\&&+\frac{1}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm1)}+2t_{k+1,2(j\pm2)}-t_{k+1,2(j\pm2)}\right)
     * \f}
     * Please note, that all values of gridpoints outside of the sparse grid are treated as 0. The following
     * picture depicts all involved grid points and temp values in order to calculate a specific point:
     * \image html prewavelets_up.png "All involved gridpoint for the up algorithm (red) and temp points between grid points (green). The gray line indicates the support of the prewavelet."
     */
    class LaplaceUpGradientPrewavelet {
      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;
        /// Pointer to sg::base::GridStorage object
        sg::base::GridStorage* storage;

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        LaplaceUpGradientPrewavelet(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        ~LaplaceUpGradientPrewavelet();

        /**
         * This operations performs the calculation of upGradient in the direction of dimension <i>dim</i>
         *
         * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result sg::base::DataVector that contains the result of the down operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        void operator()(sg::base::DataVector& source, sg::base::DataVector& result,
                        grid_iterator& index, size_t dim);
    };

  }
}

#endif /* LAPLACEUPGRADIENTMODLINEAR_HPP */
