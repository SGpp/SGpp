/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Stefanie Schraufstetter 8schraufs@in.tum.de)

#ifndef VARIABLEDISCOUNTFACTOR_HPP
#define VARIABLEDISCOUNTFACTOR_HPP

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    /**
     * This class provides a method for discounting when the risk-free rate corresponds to an axis of the grid
     */
    class VariableDiscountFactor {

      private:
        /// bounding box of the grid
        SGPP::base::BoundingBox* myBoundingBox;
        /// SGPP::base::Grid Storage object
        SGPP::base::GridStorage* storage;
        /// dimension of the risk-free rate (HW axis)
        int dim_r;

      public:
        /**
         * Std-Constructor
         *
         * @param storage the grid's storage object; needed to determine the bounding box and to iterate of the entries in the coefficient vector
         * @param dim_r dimension of the risk-free rate (HW axis)
         */
        VariableDiscountFactor(SGPP::base::GridStorage* storage, int dim_r);

        /**
         * Std-Destructor
         */
        ~VariableDiscountFactor();
        /**
         * get a vector which contains all the factors needed to multiply with another vector
         *@param factor the vector that should be calculated to multiply with another vector
         *@param T timestepsize
         */
        void getDiscountFactor(SGPP::base::DataVector& factor, double T);

    };

  }
}

#endif /* VARIABLEDISCOUNTFACTOR_HPP */
