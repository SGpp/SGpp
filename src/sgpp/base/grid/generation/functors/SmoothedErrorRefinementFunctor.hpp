/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)

#ifndef SMOOTHEDERRORREFINEMENTFUNCTOR_HPP
#define SMOOTHEDERRORREFINEMENTFUNCTOR_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/exception/application_exception.hpp"

namespace sg {
  namespace base {

    /**
     * A refinement functor, refining according to the maximal absolute values in a DataVector provided.
     * @version $HEAD$
     */
    class SmoothedErrorRefinementFunctor : public RefinementFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
         * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points)
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold
         */
        SmoothedErrorRefinementFunctor(DataVector* alpha, Grid* grid, size_t refinements_num = 1, double threshold = 0.0);

        /**
         * Destructor
         */
        virtual ~SmoothedErrorRefinementFunctor();

        virtual double operator()(GridStorage* storage, size_t seq);

        virtual double start();

        size_t getRefinementsNum();

        double getRefinementThreshold();

        void setTrainDataset(DataMatrix* trainDataset);

        void setClasses(DataVector* classes);

      protected:
        /// pointer to the vector that stores the alpha values
        DataVector* alpha;

        /// number of grid points to refine
        size_t refinements_num;

        /// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
        double threshold;

        Grid* grid;
        DataMatrix* trainDataset;
        DataVector* classes;
    };

  }
}

#endif /* SMOOTHEDERRORREFINEMENTFUNCTOR_HPP */
