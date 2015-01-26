/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef DPHIPHIDOWNBBLINEARSTRETCHEDBOUNDARY_HPP
#define DPHIPHIDOWNBBLINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/finance/basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    /**
     * Implementation of SGPP::base::sweep operator (): 1D Down for
     * Bilinearform \f$\int_{x} \frac{\partial \phi(x)}{x} \phi(x) dx\f$
     * on linear boundary grids
     */
    class DPhiPhiDownBBLinearStretchedBoundary : public DPhiPhiDownBBLinearStretched {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        DPhiPhiDownBBLinearStretchedBoundary(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~DPhiPhiDownBBLinearStretchedBoundary();

        /**
         * This operations performs the calculation of down in the direction of dimension <i>dim</i>
         *
         * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
         * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
         * result)
         * So please assure that both functions do exist!
         *
         * @param source SGPP::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result SGPP::base::DataVector that contains the result of the down operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim);
    };

    // namespace detail

  } // namespace SGPP
}

#endif /* DPHIPHIDOWNBBLINEARSTRETCHEDBOUNDARY_HPP */
