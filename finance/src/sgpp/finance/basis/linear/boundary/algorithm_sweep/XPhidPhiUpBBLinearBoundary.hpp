/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef XPHIDPHIUPBBLINEARBOUNDARY_HPP
#define XPHIDPHIUPBBLINEARBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    /**
     * Implementation of sweep operator (): 1D Up for
     * Bilinearform \f$\int_{x} x \phi(x) \frac{\partial \phi(x)}{x} dx\f$
     * on linear boundary grids
     */
    class XPhidPhiUpBBLinearBoundary : public XPhidPhiUpBBLinear {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        XPhidPhiUpBBLinearBoundary(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~XPhidPhiUpBBLinearBoundary();

        /**
         * This operations performs the calculation of up in the direction of dimension <i>dim</i>
         *
         * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
         * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
         * result)
         * So please assure that both functions do exist!
         *
         * @param source SGPP::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result SGPP::base::DataVector that contains the result of the up operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim);
    };

    // namespace detail

  } // namespace SGPP
}

#endif /* XPHIDPHIUPBBLINEARBOUNDARY_HPP */
