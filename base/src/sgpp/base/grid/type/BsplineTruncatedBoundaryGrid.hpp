// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINETRUNCATEDBOUNDARYGRID_HPP
#define BSPLINETRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Grid with Bspline basis functions with boundaries, pentagon cut
     * @todo (pflueged) include for factory exception is missing in several classes which use it. It only works, as it is include by a header loaded previously.
     */
    class BsplineTruncatedBoundaryGrid : public Grid {
      protected:
        BsplineTruncatedBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified bspline basis functions
         *
         * @param dim the dimension of the grid
         * @param degree the bspline's degree
         */
        BsplineTruncatedBoundaryGrid(size_t dim, size_t degree);

        /**
         * Destructor
         */
        virtual ~BsplineTruncatedBoundaryGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

        virtual void serialize(std::ostream& ostr);
        virtual size_t getDegree();

      protected:
        // degree of Bspline
        size_t degree;

        const SBsplineBoundaryBase* basis_;


    };

  }
}

#endif /* BSPLINETRUNCATEDBOUNDARYGRID_HPP */
