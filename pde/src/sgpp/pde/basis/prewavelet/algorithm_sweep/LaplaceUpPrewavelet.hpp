/* ****************************************************************************
* Copyright (C) 2008-2009 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de)
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LAPLACEUPPREWAVELET_HPP
#define LAPLACEUPPREWAVELET_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sg {
  namespace pde {



    /**
     * Implements the down Method needed for the Laplace operator on prewavelet grids.
     * Here, the prewavlets have an advantage over the normal linear base due to the semi
     * othogonality. That means, for the calculation, only the gridpoint in the same level have
     * to be taken into account. The matrix entries are the following:
     * \f{eqnarray*}{
     * m_{(1,1)(1,1)}&=&\frac{1}{3}\\m_{(2,1)(2,3)}&=&\frac{1}{25}\\m_{(l,1)(l,1)}&=&m_{(l,2^{l}-1)(l,2^{l}-1)}=\frac{44}{75}h_{l}\\m_{(l,1)(l,3)}&=&m_{(l,2^{l}-1)(l,2^{l}-3)}=\frac{11}{75}h_{l}\\m_{(l,i)(l,i)}&=&\frac{18}{25}h_{l}\\m_{(l,i)(l,i\pm2)}&=&\frac{2}{15}h_{l}\\m_{(l,i)(l,i\pm4)}&=&-\frac{11}{75}h_{l}
     * \f}
     * With that, the calculation of the entire matrix is completed, thus no additional up-part is needed.
     */
    class LaplaceUpPrewavelet {
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
        LaplaceUpPrewavelet(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        ~LaplaceUpPrewavelet();

        /**
         * This operations performs the calculation of down in the direction of dimension <i>dim</i>
         *
         * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result sg::base::DataVector that contains the result of the up operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        void operator()(sg::base::DataVector& source, sg::base::DataVector& result,
                        grid_iterator& index, size_t dim);

      protected:

    };



  }
}

#endif /* LAPLACEUPPREWAVELET_HPP */
