/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)

#ifndef WEIGHTEDERRORREFINEMENTFUNCTOR_HPP
#define WEIGHTEDERRORREFINEMENTFUNCTOR_HPP

#include "sgpp_base.hpp"
#include <limits>
#include "base/exception/application_exception.hpp"

namespace sg {
  namespace base {

    /**
     * A refinement functor, refining according to the maximal absolute values in a DataVector provided.
     * @version $HEAD$
     */
    class WeightedErrorRefinementFunctor : public RefinementFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
         * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points)
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold
         */
        WeightedErrorRefinementFunctor(DataVector* alpha, Grid* grid, size_t refinements_num = 1, double threshold = 0.0);

        /**
         * Destructor
         */
        virtual ~WeightedErrorRefinementFunctor();


        /*
         * Calculates
         * s_j = \sum_{i=0}^{n} phi_j(x_i) * alpha_j * err_i^2
         *
         * where
         * j = index of the current basis function
         * s_j = indicator value
         *
         * u(x) = interpolant
         *
         * i = index of the current data point
         * x_i = x value of data point
         * y_i = y value of data point
         * err_i = u(x_i) - y_i
         *
         * n = number of data points
         */
        virtual double operator()(GridStorage* storage, size_t seq);

        virtual double start();

        size_t getRefinementsNum();

        double getRefinementThreshold();

        void setTrainDataset(DataMatrix* trainDataset);

        void setClasses(DataVector* classes);

        void setErrors(DataVector* classes);

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
    };

  }
}

#endif /* WEIGHTEDERRORREFINEMENTFUNCTOR_HPP */
