// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODFUNDAMENTALSPLINEGRID_HPP
#define MODFUNDAMENTALSPLINEGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Grid with modified fundamental spline basis functions
     */
    class ModFundamentalSplineGrid : public Grid {
      protected:
        /**
         * This constructor creates a new GridStorage out of the stream.
         *
         * @param istr inputstream that contains the grid information
         */
        ModFundamentalSplineGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified fundamental spline basis functions
         *
         * @param dim the dimension of the grid
         * @param degree fundamental spline degree
         */
        ModFundamentalSplineGrid(size_t dim, size_t degree);

        /**
         * Destructor.
         */
        virtual ~ModFundamentalSplineGrid();

        /**
         * @return string that identifies the grid type uniquely
         */
        virtual SGPP::base::GridType getType();

        /**
         * @return fundamental spline basis
         */
        virtual const SBasis& getBasis();

        /**
         * @return pointer to a GridGenerator object
         */
        virtual GridGenerator* createGridGenerator();

        /**
         * reads a grid out of a string
         *
         * @param istr string that contains the grid information
         * @return grid
         */
        static Grid* unserialize(std::istream& istr);

        /**
         * Serializes the grid.
         *
         * @param ostr stream to which the grid is written
         */
        virtual void serialize(std::ostream& ostr);

        /**
         * @return B-spline degree
         */
        virtual size_t getDegree();

      protected:
        /// fundamental spline degree
        size_t degree;
        /// fundamental spline basis
        const SFundamentalSplineModifiedBase* basis_;
    };

  }
}

#endif /* MODFUNDAMENTALSPLINEGRID_HPP */
