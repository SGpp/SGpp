/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef DIRICHLETUPDATEVECTOR_HPP
#define DIRICHLETUPDATEVECTOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sg {
  namespace base {

    /**
     * This class is useful if you do some PDE calculations with Dirichlet Boundary
     * Conditions. Doing this, e.g. you might wish to add some solution from a timestep to
     * the current coefficients of the sparse grid. If you are using Dirichlet
     * conditions you mustn't overwrite the values on the boundaries in your coefficient
     * vector.
     *
     * This class implements a functor that uses the Bounding Box of the grid to determine, if
     * a boundary has to implement Dirichlet boundary conditions. In that case, simply use this
     * to replace all values in the update vector on these boundaries with zero, so you can
     * safely add the resulting vector to your solution.
     */
    class DirichletUpdateVector {
      private:
        //  /// bounding box of the grid
        //  BoundingBox* myBoundingBox;
        //  /// stretching of the grid
        //  Stretching* myStretching;
        /// Grid Storage object
        GridStorage* storage;

      public:
        /**
         * Std-Constructor
         *
         * @param storage the grid's storage object; needed to determine the bounding box and to iterate of the entries in the coefficient vector
         */
        DirichletUpdateVector(GridStorage* storage);

        /**
         * Std-Destructor
         */
        ~DirichletUpdateVector();

        /**
         * Replace the boundary entries in updateVector with the one from sourceVector
         * only in that dimension, for which Dirichlet Boundary Conditions
         * were specified
         *
         * @param updateVector the vector that should be updated
         * @param sourceVector the vector that contains the correct boundary values
         */
        void applyDirichletConditions(DataVector& updateVector, DataVector& sourceVector);

        /**
         * Replace the boundary entries in updateVector with Zero only in that dimension, for which Dirichlet Boundary Conditions
         * were specified
         *
         * @param updateVector the vector that should be updated
         */
        void setBoundariesToZero(DataVector& updateVector);

        /**
         * Replace the inner entries in updateVector with Zero only in that dimension, for which Dirichlet Boundary Conditions
         * were specified
         *
         * @param updateVector the vector that should be updated
         */
        void setInnerPointsToZero(DataVector& updateVector);
        /**
         * Multiplies the values on the boundary with vector
         * @param updateVector the vector that should be updated
         * @param factor the vector contains corresponding values
         */
        void multiplyBoundaryVector(DataVector& updateVector, DataVector& factor);
        /**
        * Multiplies the values on the boundary with a constant value
        *
        * @param updateVector the vector that should be updated
        * @param value the value that is multiplied with the value on the boundaries
        */
        void multiplyBoundary(DataVector& updateVector, double value);

        /**
         * Multiplies the values of the points in the vector that meet the predicate condition by the constant value.
         * Calling this method with a function pointer that returns true if point->isInnerPoint() and false otherwise gives the same result as
         * the multiplyBoundary method.
         */
        void multiply(DataVector& updateVector, double value, bool (*predicate)(GridIndex*, GridStorage*));
    };

  }
}

#endif /* DIRICHLETUPDATEVECTOR_HPP */
