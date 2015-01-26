/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef OPERATIONLTWODOTPRODUCTLINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONLTWODOTPRODUCTLINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/pde/algorithm/StdUpDown.hpp>


namespace sg {
  namespace pde {

    /**
     * Implements the standard L 2 scalar product on linear boundary grids
     *
     * @version $HEAD$
     */
    class OperationLTwoDotProductLinearStretchedBoundary: public StdUpDown {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationLTwoDotProductLinearStretchedBoundary(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLTwoDotProductLinearStretchedBoundary();

      protected:
        /**
         * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the up-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        /**
         * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
         * Applies the down-part of the one-dimensional mass matrix in one dimension.
         * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLTWODOTPRODUCTLINEARSTRETCHEDBOUNDARY_HPP */
