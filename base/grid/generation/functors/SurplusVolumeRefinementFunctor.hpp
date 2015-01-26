/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef SURPLUSVOLUMEREFINEMENTFUNCTOR_HPP
#define SURPLUSVOLUMEREFINEMENTFUNCTOR_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * A refinement functor, refining according to the maximal absolute values in a DataVector provided,
     * weighted with the corresponding basis function's surplus, i.e., with @f$2^{-|\vec{l}|_1} = 2^{\sum_{k=1}^d l_d}@f$.
     * @version $HEAD$
     * @todo (pflueged) add unit tests
     */
    class SurplusVolumeRefinementFunctor : public RefinementFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
         * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points), default: 1
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold, default: 0.0
         */
        SurplusVolumeRefinementFunctor(DataVector* alpha, size_t refinements_num = 1, double threshold = 0.0);

        /**
         * Destructor
         */
        virtual ~SurplusVolumeRefinementFunctor();

        virtual double operator()(GridStorage* storage, size_t seq);

        virtual double start();

        size_t getRefinementsNum();

        double getRefinementThreshold();

      protected:
        /// pointer to the vector that stores the alpha values
        DataVector* alpha;

        /// number of grid points to refine
        size_t refinements_num;

        /// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
        double threshold;
    };

  }
}

#endif /* SURPLUSVOLUMEREFINEMENTFUNCTOR_HPP */
