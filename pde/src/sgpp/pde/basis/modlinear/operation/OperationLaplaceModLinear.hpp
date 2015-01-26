/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACEMODLINEAR_HPP
#define OPERATIONLAPLACEMODLINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

namespace sg {
  namespace pde {

    /**
     * Implementation of Laplace for mod linear functions
     *
     * @version $HEAD$
     */
    class OperationLaplaceModLinear : public UpDownOneOpDim {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationLaplaceModLinear(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceModLinear();

      protected:
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACEMODLINEAR_HPP */
