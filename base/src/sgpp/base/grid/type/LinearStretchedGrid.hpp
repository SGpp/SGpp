/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pfueged@in.tum.de)

#ifndef LINEARSTRETCHEDGRID_HPP
#define LINEARSTRETCHEDGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linearstretched base functions
     */
    class LinearStretchedGrid : public Grid {
      protected:
        LinearStretchedGrid(std::istream& istr);

      public:
        /**
         * Constructor LinearStretched Grid without boundaries
         *
         * @param dim the dimension of the grid
         */
        LinearStretchedGrid(size_t dim);

        /**
         * Constructor LinearStretched Grid
         *
         * @param BB the BoundingBox of the grid
         */
        LinearStretchedGrid(Stretching& BB);

        /**
         * Destructor
         */
        virtual ~LinearStretchedGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();
        static Grid* unserialize(std::istream& istr);
    };

  }
}

#endif /* LINEARSTRETCHEDGRID_HPP */
