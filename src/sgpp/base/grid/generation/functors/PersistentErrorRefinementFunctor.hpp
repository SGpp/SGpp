/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)

#ifndef PERSISTENTERRORREFINEMENTFUNCTOR_HPP
#define PERSISTENTERRORREFINEMENTFUNCTOR_HPP

#include <vector>
#include "sgpp_base.hpp"

namespace sg {
  namespace base {

    /**
     * A refinement functor, refining according to the maximal absolute values in a DataVector provided.
     * @version $HEAD$
     */
    class PersistentErrorRefinementFunctor : public RefinementFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
         * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points)
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold
         */
        PersistentErrorRefinementFunctor(DataVector* alpha, Grid* grid, size_t refinements_num = 1, double threshold = 0.0);

        /**
         * Destructor
         */
        virtual ~PersistentErrorRefinementFunctor();


        /*
         * Uses a discounting error indicator vector to retrieve the value.
         *
         * The error indicator is updated using the rule
         * disc_err_{n+1} = disc_err_{n} * BETA + weight_error_per_basis * (1-BETA)
         * where 0 < BETA < 1
         */
        virtual double operator()(GridStorage* storage, size_t seq);

        virtual double start();

        size_t getRefinementsNum();

        double getRefinementThreshold();

        void setTrainDataset(DataMatrix* trainDataset_);
        void setClasses(DataVector* classes_);
        void setErrors(DataVector* errors);

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
        DataVector* errors;
        DataVector* accum;
    };

  }
}

#endif /* PERSISTENTERRORREFINEMENTFUNCTOR_HPP */
