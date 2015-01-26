/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SURPLUSREFINEMENTFUNCTOR_HPP
#define SURPLUSREFINEMENTFUNCTOR_HPP

#include <iostream>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * A refinement functor, refining according to the maximal absolute values in a DataVector provided.
     * @version $HEAD$
     */
    class SurplusRefinementFunctor : public RefinementFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
         * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points)
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold
         */
        SurplusRefinementFunctor(DataVector* alpha, size_t refinements_num = 1, double threshold = 0.0);

        /**
         * Destructor
         */
        virtual ~SurplusRefinementFunctor();

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

#endif /* SURPLUSREFINEMENTFUNCTOR_HPP */
