/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef COARSENINGFUNCTOR_HPP
#define COARSENINGFUNCTOR_HPP

#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * Abstract class that defines the interfaces that coarsening functors have to provide.
     * @version $HEAD$
     */
    class CoarseningFunctor {
      public:
        typedef double value_type;

        /**
         * Constructor
         */
        CoarseningFunctor() {}

        /**
         * Destructor
         */
        virtual ~CoarseningFunctor() {}

        /**
         * This should be returning a coarsening value for every grid point.
         * The point with the lowest value will be removed first.
         *
         * @param storage pointer to the grids storage object
         * @param seq sequence number in the coefficients array
         *
         * @return refinement value
         */
        virtual double operator()(GridStorage* storage, size_t seq) = 0;

        /**
         * This should return the initial value of coarsening criterion (e.g. alpha or error).
         *
         * @return the initial value
         */
        virtual double start() = 0;

        /**
         * Returns the maximal number of points that should be removed.
         *
         * The maximal number of points to removed is set in the constructor of implementation class.
         *
         * @return number of points that should removed. Default value: 1.
         */
        virtual size_t getRemovementsNum() {
          return 1;
        }

        /**
         * Returns the threshold value.
         *
         * Only the grid points with absolute value of coarsening criterion (e.g. alpha) less
         * or equal to this threshold will be removed
         *
         * @return threshold value for refinement. Default value: 0.
         */
        virtual double getCoarseningThreshold() = 0;
    };

  }
}

#endif /* COARSENINGFUNCTOR_HPP */
