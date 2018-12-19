// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#ifndef SGPP_GRIDPOINTBASEDCOARSENINGFUNCTOR_HPP
#define SGPP_GRIDPOINTBASEDCOARSENINGFUNCTOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/functors/MultiGridCoarseningFunctor.hpp>

#include <vector>
#include <map>
#include <string>


namespace sgpp {
    namespace datadriven {

/**
 * Grid Point based coarsening for classification problems solved by
 * a SG density estimation approach. The scoring is only based on
 * function values of the respective PDF at grid points.
 */
        class GridPointBasedCoarseningFunctor : public MultiGridCoarseningFunctor {
        public:
            /**
             * Constructor.
             *
             * @param grids Vector of grids. current_grid_index specifies the grid to be coarsened
             * @param alphas Vector of surpluses related to the grids
             * @param coarsenings_num Maximum number of coarsenings done
             * @param pre_compute Flag for precomputation of necessary grid evals. If true preComputeEvaluations needs to be called before each coarsening step
             * @param threshold Threshold for coarsening scores
             */
            GridPointBasedCoarseningFunctor(std::vector<base::Grid*> grids,
                                            std::vector<base::DataVector*> alphas,
                                            size_t coarsenings_num = 1,
                                            bool pre_compute = false,
                                            double threshold = 0.0);

            double operator()(base::GridStorage& storage,
                              size_t seq) const override;
            double start() const override;
            size_t getRemovementsNum() const override;
            double getCoarseningThreshold() const override;
            virtual ~GridPointBasedCoarseningFunctor() {}

            void setGridIndex(size_t grid_index) override;
            size_t getNumGrids() override;

            void preComputeEvaluations() override;


        protected:
            std::vector<base::Grid*> grids;
            std::vector<base::DataVector*> alphas;

            size_t current_grid_index;
            size_t coarsenings_num;
            double threshold;

            bool pre_compute;

            /**
             * Stores grid evaluations at all grids (vector) at the union
             * of grid point coordinates over all grids (hashed by string representation
             * in the map)
             */
            std::vector<std::map<std::string, double>> pre_comp_evals;
        };
    }  // namespace datadriven
}  // namespace sgpp

#endif /* GRIDPOINTBASEDCOARSENINGFUNCTOR_HPP */
