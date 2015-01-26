/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef ANOVACOARSENINGFUNCTOR_HPP_
#define ANOVACOARSENINGFUNCTOR_HPP_

#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>


namespace sg {
  namespace base {

    /*
     * Coarsening functor that implements the ANOVA compression scheme. It calculates
     * the total Variance as a sum of L2 Norms of all ANOVA components and keeps only
     * those functions that lie in the ANOVA components with top [thershold]% of the
     * total variance
     */
    class ANOVACoarseningFunctor : public CoarseningFunctor {
      public:
        /**
         * Constructor.
         *
         * @param alpha DataVector that is basis for coarsening decisions. The i-th entry corresponds to the i-th grid point.
         * @param removements_num Number of grid points which should be removed (if possible - there could be less removable grid points)
         * @param threshold The absolute value of the entries have to be greater or equal than the threshold
         * @param storage grid storage
         */
        ANOVACoarseningFunctor(DataVector* alpha, size_t removements_num, double threshold, GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~ANOVACoarseningFunctor();


        /**
         * This should be returning a coarsening value for every grid point.
         * The point with the lowest value will be removed first.
         *
         * @param storage pointer to the grids storage object
         * @param seq sequence number in the coefficients array
         *
         * @return refinement value
         */
        virtual double operator()(GridStorage* storage, size_t seq);

        /**
         * This should return the initial value of coarsening criterion (e.g. alpha or error).
         *
         * @return the initial value
         */
        virtual double start();


        /**
         * Returns the maximal number of points that should be removed.
         *
         * The maximal number of points to removed is set in the constructor of implementation class.
         *
         * @return number of points that should removed. Default value: 1.
         */
        size_t getRemovementsNum();


        /**
         * Returns the threshold value.
         *
         * Only the grid points with absolute value of coarsening criterion (e.g. alpha) less
         * or equal to this threshold will be removed
         *
         * @return threshold value for refinement. Default value: 0.
         */
        double getCoarseningThreshold();


      protected:

        /**
         * Returns the index of the ANOVA component the basis function with the
         * give index lies in
         *
         * @param index index of the basis functino
         * @return anova component index
         */
        int getANOVAComponentIndex(GridStorage::index_type& index);

        /// pointer to the vector that stores the alpha values
        DataVector* alpha;

        /// number of grid points to remove
        int removements_num;

        /// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
        double threshold;





      private:

        typedef struct ANOVA_Values {
          int component_index;
          double value;
          double refinement_value;
          int points_num;

        } tANOVAValues;

        struct Sorter {
          bool operator() (tANOVAValues i, tANOVAValues j) {
            return (i.value > j.value);
          }
        }
        sorter;

        std::vector<tANOVAValues> anova_variances;
        tANOVAValues** anova_variances_pointers;
    };

  } /* namespace base */
} /* namespace sg */
#endif /* ANOVACOARSENINGFUNCTOR_HPP_ */
