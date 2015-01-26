/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef DIRICHLETGRIDCONVERTER_HPP
#define DIRICHLETGRIDCONVERTER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

namespace sg {
  namespace base {

    /**
     * This class handles the conversion of a boundary grid with dirichlet
     * boundary conditions into an inner. This is useful in case of solving a linear
     * system of the grid's ansatzfunctions' coefficients: Here the boundary points
     * are needed during solving the system (they constant because dirichlet boundary
     * conditions), so a lot of work can be saved if only the inner points are used.
     *
     * @version $HEAD$
     */
    class DirichletGridConverter {
      private:
        /// number of the boundary grid's grid points
        size_t numTotalGridPoints;
        /// number of the inner grid's grid points
        size_t numInnerGridPoints;
        /// array to store the position of i-th inner point in the boundary grid's coefficients
        size_t* conCoefArray;
        /// ensure that buildInnerGridWithCoefs can only be called once
        bool bFirstTime;

      public:

        /**
         * Constructor
         */
        DirichletGridConverter();

        /**
         * Destructor
         */
        ~DirichletGridConverter();

        /**
         * builds a sparse grid without any boundaries from a sparse grid with boundaries. In addition
         * the coefficient vector is also created and initialized for the inner grid.
         *
         * @param BoundaryGrid the boundary grid whose inner should be extracted
         * @param BoundaryCoefs the boundary's grid coefficients
         * @param InnerGrid Pointer to the inner grid, initialized in this method
         * @param InnerCoefs Pointer to the inner grid's coefficients, initialized in this mehtod
         */
        void buildInnerGridWithCoefs(Grid& BoundaryGrid, DataVector& BoundaryCoefs, Grid** InnerGrid, DataVector** InnerCoefs);

        /**
         * rebuilds a sparse grid without any boundaries from a sparse grid with boundaries. In addition
         * the coefficient vector is also created and initialized for the inner grid.
         *
         * @param BoundaryGrid the boundary grid whose inner should be extracted
         * @param BoundaryCoefs the boundary's grid coefficients
         * @param InnerGrid Pointer to the inner grid, initialized in this method
         * @param InnerCoefs Pointer to the inner grid's coefficients, initialized in this mehtod
         */
        void rebuildInnerGridWithCoefs(Grid& BoundaryGrid, DataVector& BoundaryCoefs, Grid** InnerGrid, DataVector** InnerCoefs);

        /**
         * copies the inner grid's coefficients to the identical (inner) ansatzfunctions in the boundary grid.
         *
         * Here a very check, due to performance issues, is implemented to sync
         * both vectors by checking the size of both vectors. It must
         * match to the creation size determined in buildInnerGridWithCoefs.
         *
         * @param BoundaryCoefs the boundary grid's coefficients
         * @param InnerCoefs the inner grid's coefficients
         */
        void updateBoundaryCoefs(DataVector& BoundaryCoefs, DataVector& InnerCoefs);

        /**
         * copies the boundary grid's inner coefficients to the inner grid coefficients
         *
         * Here a very check, due to performance issues, is implemented to sync
         * both vectors by checking the size of both vectors. It must
         * match to the creation size determined in buildInnerGridWithCoefs.
         *
         * @param BoundaryCoefs the boundary grid's coefficients
         * @param InnerCoefs the inner grid's coefficients
         */
        void calcInnerCoefs(DataVector& BoundaryCoefs, DataVector& InnerCoefs);
    };

  }
}

#endif /* DIRICHLETGRIDCONVERTER_HPP */
