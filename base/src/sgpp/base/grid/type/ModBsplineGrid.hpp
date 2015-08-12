// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODBSPLINEGRID_HPP
#define MODBSPLINEGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Grid with modified Bspline basis functions
     */
    class ModBsplineGrid : public Grid {
      protected:
        /**
         * This constructor creates a new GridStorage out of the stream.
         *
         * @param istr inputstream that contains the grid information
         */
        ModBsplineGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified bspline basis functions
         *
         * @param dim the dimension of the grid
           * @param degree the bspline's degree
         */
        ModBsplineGrid(size_t dim, size_t degree);

        /**
         * Destructor.
         */
        virtual ~ModBsplineGrid();

        /**
         * @return string that identifies the grid type uniquely
         */
        virtual const char* getType();

        /**
         * @return B-spline basis
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
        /// B-spline degree
        size_t degree;
        /// B-spline basis
        const SBsplineModifiedBase* basis_;


    };

  }
}

#endif /* MODBSPLINEGRID_HPP */
