// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PHIPHIUPBBLINEARSTRETCHEDBOUNDARY_HPP
#define PHIPHIUPBBLINEARSTRECHEDBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {



    /**
     * Implementation of SGPP::base::sweep operator (): 1D Up for
     * Bilinearform \f$\int_{x} \phi(x) \phi(x) dx\f$
     * on linear boundary grids
     */
    class PhiPhiUpBBLinearStretchedBoundary : public PhiPhiUpBBLinearStretched {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        PhiPhiUpBBLinearStretchedBoundary(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~PhiPhiUpBBLinearStretchedBoundary();

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

#endif /* PHIPHIUPBBLINEARSTRETCHEDBOUNDARY_HPP */