// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef XPHIDPHIUPBBLINEARSTRETCHED_HPP
#define XPHIDPHIUPBBLINEARSTRETCHED_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {



    /**
     * Implementation of sweep operator (): 1D Up for
     * Bilinearform \f$\int_{x} x \phi(x) \frac{\partial \phi(x)}{x} dx\f$
     */
    class XPhidPhiUpBBLinearStretched {
      protected:
        typedef SGPP::base::GridStorage::grid_iterator grid_iterator;

        /// Pointer to SGPP::base::GridStorage object
        SGPP::base::GridStorage* storage;
        /// Pointer to the stretching Obejct
        SGPP::base::Stretching* stretching;

      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        XPhidPhiUpBBLinearStretched(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        ~XPhidPhiUpBBLinearStretched();

        /**
         * This operations performs the calculation of up in the direction of dimension <i>dim</i>
         * on a a grid with fix Dirichlet 0 boundary condition
         *
         * @param source SGPP::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result SGPP::base::DataVector that contains the result of the up operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        void operator()(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim);

      protected:

        /**
         * recursive function for the calculation of Up
         *
         * @param source SGPP::base::DataVector that contains the coefficients of the ansatzfunction
         * @param result SGPP::base::DataVector in which the result of the operation is stored
         * @param index reference to a griditerator object that is used navigate through the grid
         * @param dim the dimension in which the operation is executed
         * @param fl function value on the left boundary, reference parameter
         * @param fr function value on the right boundary, reference parameter
         */
        void rec(SGPP::base::DataVector& source, SGPP::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr);

    };

    // namespace detail

  } // namespace SGPP
}

#endif /* XPHIDPHIUPBBLINEARSTRETCHEDBOUNDARY_HPP */