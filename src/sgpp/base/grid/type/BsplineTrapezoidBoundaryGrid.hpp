/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_GRID_TYPE_BSPLINETRAPEZOIDBOUNDARYGRID_HPP
#define SGPP_BASE_GRID_TYPE_BSPLINETRAPEZOIDBOUNDARYGRID_HPP

#include <iostream>

#include "base/grid/Grid.hpp"

namespace sg {
  namespace base {

    /**
     * Boundary grid with B-spline basis functions.
     */
    class BsplineTrapezoidBoundaryGrid : public Grid {
      public:
        /**
         * Constructor.
         *
         * @param dim       number of dimensions
         * @param degree    B-spline degree
         */
        BsplineTrapezoidBoundaryGrid(size_t dim, size_t degree);

        /**
         * Destructor.
         */
        virtual ~BsplineTrapezoidBoundaryGrid();

        /**
         * @return  identifying grid type string
         */
        virtual const char* getType();

        /**
         * @return grid generator for this grid type
         */
        virtual GridGenerator* createGridGenerator();

        /**
         * @param[out] ostr     output stream as target of serialization
         */
        virtual void serialize(std::ostream& ostr);

        /**
         * @param istr  input stream containing the serialization
         * @return      pointer to newly generated deserialized grid
         */
        static Grid* unserialize(std::istream& istr);

        /**
         * @return  B-spline degree
         */
        virtual size_t getDegree();

      protected:
        /// B-spline degree
        size_t degree;

        /**
         * Deserialization constructor.
         *
         * @param istr  serialized grid
         */
        BsplineTrapezoidBoundaryGrid(std::istream& istr);
    };

  }
}

#endif
