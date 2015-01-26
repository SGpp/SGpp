/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACEMODLINEAR_HPP
#define OPERATIONLAPLACEMODLINEAR_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
         * @param storage the grid's SGPP::base::GridStorage object
         */
        OperationLaplaceModLinear(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceModLinear();

      protected:
        virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACEMODLINEAR_HPP */
