// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PERSISTENTERRORREFINEMENTFUNCTOR_HPP
#define PERSISTENTERRORREFINEMENTFUNCTOR_HPP

#include <vector>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * A refinement functor, refining according to the maximal absolute values in a DataVector provided.
     */
    class PersistentErrorRefinementFunctor : public RefinementFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for refinement decisions. The i-th entry corresponds to the i-th grid point.
         * @param grid grid
         * @param refinements_num Number of grid points which should be refined (if possible - there could be less refinable grid points)
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold
         */
        PersistentErrorRefinementFunctor(DataVector* alpha, Grid* grid, size_t refinements_num = 1, float_t threshold = 0.0);

        /**
         * Destructor
         */
        virtual ~PersistentErrorRefinementFunctor() override;


        /*
         * Uses a discounting error indicator vector to retrieve the value.
         *
         * The error indicator is updated using the rule
         * disc_err_{n+1} = disc_err_{n} * BETA + weight_error_per_basis * (1-BETA)
         * where 0 < BETA < 1
         *
         * @param storage pointer to the grids storage object
         * @param seq sequence number in the coefficients array
         *
         * @return refinement value
         */
        virtual float_t operator()(GridStorage* storage, size_t seq) override;

        virtual float_t start() override;

        size_t getRefinementsNum();

        float_t getRefinementThreshold();

        void setTrainDataset(DataMatrix* trainDataset_);
        void setClasses(DataVector* classes_);
        void setErrors(DataVector* errors);

      protected:
        /// pointer to the vector that stores the alpha values
        DataVector* alpha;

        /// number of grid points to refine
        size_t refinements_num;

        /// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
        float_t threshold;

        Grid* grid;
        DataMatrix* trainDataset;
        DataVector* classes;
        DataVector* errors;
        DataVector* accum;
    };

  }
}

#endif /* PERSISTENTERRORREFINEMENTFUNCTOR_HPP */
