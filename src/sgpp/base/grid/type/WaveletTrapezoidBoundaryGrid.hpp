/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_GRID_TYPE_WAVELETTRAPEZOIDBOUNDARYGRID_HPP
#define SGPP_BASE_GRID_TYPE_WAVELETTRAPEZOIDBOUNDARYGRID_HPP

#include <iostream>

#include "base/grid/Grid.hpp"

namespace sg {
  namespace base {

    /**
     * Boundary grid with wavelet basis functions.
     */
    class WaveletTrapezoidBoundaryGrid : public Grid {
      public:
        /**
         * Constructor.
         *
         * @param dim       number of dimensions
         */
        WaveletTrapezoidBoundaryGrid(size_t dim);

        /**
         * Destructor.
         */
        virtual ~WaveletTrapezoidBoundaryGrid();

        /**
         * @return  identifying grid type string
         */
        virtual const char* getType();

        /**
         * @return grid generator for this grid type
         */
        virtual GridGenerator* createGridGenerator();

        /**
         * @param istr  input stream containing the serialization
         * @return      pointer to newly generated deserialized grid
         */
        static Grid* unserialize(std::istream& istr);

      protected:
        /**
         * Deserialization constructor.
         *
         * @param istr  serialized grid
         */
        WaveletTrapezoidBoundaryGrid(std::istream& istr);
    };

  }
}

#endif
