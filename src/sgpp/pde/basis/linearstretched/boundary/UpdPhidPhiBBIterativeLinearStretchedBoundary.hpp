/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef UPDPHIDPHIDOWNBBITERATIVELINEARSTRETCHEDBOUNDARY_HPP
#define UPDPHIDPHIDOWNBBITERATIVELINEARSTRETCHEDBOUNDARY_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {

    /**
     * This class is helper class to implement the complete Up
     * of following bilinearform \f$\int_{x} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x} dx\f$
     * for a given dimension by an iterative algorithms on adaptive
     * Sparse Grids with linear ansatzfunctions with boundaries.
     *
     * This is possible due to the fact that the operator's
     * matrix has only entries on the diagonal.
     *
     * -> the Up/Down can be implemented by iterating over
     * all ansatzfunctions
     */
    class UpdPhidPhiBBIterativeLinearStretchedBoundary {
      private:
        /// Pointer to the grid's storage object
        sg::base::GridStorage* storage;

      public:
        /**
         * Constructor
         *
         * @param storage Pointer to the grid's storage object
         */
        UpdPhidPhiBBIterativeLinearStretchedBoundary(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        ~UpdPhidPhiBBIterativeLinearStretchedBoundary();

        /**
         * This operations performs the calculation of Up in the direction of dimension <i>dim</i>
         * of following bilinearform: \f$\int_{x} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x} dx\f$
         *
         * @param alpha sg::base::DataVector that contains the gridpoint's coefficients
         * @param result sg::base::DataVector that contains the result of the down operation
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* UPDPHIDPHIDOWNBBITERATIVELINEARSTRETCHEDBOUNDARY_HPP */
