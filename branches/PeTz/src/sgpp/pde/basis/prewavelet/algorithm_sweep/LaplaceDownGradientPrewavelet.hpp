/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de)
// @author xander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LAPLACEDOWNGRADIENTPREWAVELET_HPP
#define LAPLACEDOWNGRADIENTPREWAVELET_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {



    /**
     * Implements the downGradient Method needed for the Laplace operator on prewavelet grids.
     * The calculation is done iterative and utilizes the following temp variables:
     * \f{eqnarray*}{
     * t_{k,j}&=&\frac{1}{h_{k}}\left(\frac{32}{10}u_{k,j}+\frac{8}{10}u_{k,j\pm2}\right) \quad(k,j)\in G_{n}^{1}\\ t_{k,j}&=&\frac{1}{h_{k}}\left(\frac{24}{10}u_{k,j}+\frac{8}{10}u_{k,j\pm2}\right) \quad(k,j)\in G_{n}^{1}\mbox{ and next to border}\\ t_{k,j}&=&t_{k-1,\frac{j}{2}}-\frac{1}{h_{k}}\left(\frac{23}{10}u_{k,j\pm1}+\frac{1}{10}u_{k,j\pm3}\right) \quad(k,j)\notin G_{n}^{1}\\ t_{k,j}&=&t_{k-1,\frac{j}{2}}-\frac{1}{h_{k}}\left(\frac{22}{10}u_{k,j\pm1}+\frac{1}{10}u_{k,j\pm3}\right) \quad(k,j)\notin G_{n}^{1}\mbox{ and \ensuremath{(k,j\pm1)}next to border}
     * \f}
     *
     * The calculation of the results is as follows:
     * \f{eqnarray*}{
     r_{k,j}&=&-\frac{6}{10}t_{k-1,\frac{j\pm1}{2}}\\
     &&+\frac{1}{h_{k}}u_{k,j}\cdot\left[\frac{612}{100}\mbox{ if \ensuremath{(k,j)} inner point, or }
     \frac{356}{100}\mbox{ if \ensuremath{(k,j)} border point}\right]\\
     &&+\frac{1}{h_{k}}u_{k,j\pm2}\cdot\left[\frac{256}{100} \mbox{ if \ensuremath{(k,j)} and \ensuremath{(k,j\pm2)} inner point, or }
     \frac{242}{100} \mbox{ if \ensuremath{(k,j)} or \ensuremath{(k,j\pm2)} border point, or }
     \frac{228}{100} \mbox{ if \ensuremath{(k,j)} border point}\right]\\
     &&+\frac{1}{h_{k}}u_{k,j\pm4}\cdot\frac{14}{100}\\r_{1,1}&=&4u_{1,1}
     * \f}
     * Please note, that all values of gridpoints outside of the sparse grid are treated as 0. The following
     * picture depicts all involved grid points and temp values in order to calculate a specific point:
     * \image html prewavelets_down.png "All involved gridpoint for the up algorithm (red) and temp points between grid points (green). The gray line indicates the support of the prewavelet."
     */

    class LaplaceDownGradientPrewavelet {
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
        LaplaceDownGradientPrewavelet(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        ~LaplaceDownGradientPrewavelet();

        /**
         * This operations performs the calculation of downGradient in the direction of dimension <i>dim</i>
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

#endif /* LAPLACEDOWNGRADIENTPREWAVELET_HPP */
