// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#ifndef GRIDPOINTBASEDREFINEMENTFUNCTOR_HPP
#define GRIDPOINTBASEDREFINEMENTFUNCTOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>

#include <vector>
#include <map>


namespace sgpp {
  namespace datadriven {

    /**
     * Grid Point based refinement for classification problems solved by
     * a SG density estimation approach. The scoring is only based on
     * function values of the respective PDF at grid points.
     */
    class GridPointBasedRefinementFunctor : public MultiGridRefinementFunctor {

    public:
      /**
       * Constructor.
       *
       * @param grids Vector of grids. current_grid_index specifies the grid to be refined
       * @param alphas Vector of surpluses related to the grids
       * @param refinements_num Maximum number of refinements done
       * @param level_penalize If a level penalizing is multiplied to the score (2^{|l|_1})
       * @param pre_compute Flag for precomputation of necessary grid evals. If true preComputeEvaluations needs to be called before each refinement step
       * @param threshold Threshold for refinement scores
       */
      GridPointBasedRefinementFunctor(std::vector<base::Grid*> grids,
				      std::vector<base::DataVector*> alphas,
				      size_t refinements_num = 1,
				      bool level_penalize = false,
				      bool pre_compute = false,
				      double threshold = 0.0);

      double operator()(base::GridStorage& storage,
			size_t seq) const override;
      double start() const override;
      size_t getRefinementsNum() const override;
      double getRefinementThreshold() const override;
      virtual ~GridPointBasedRefinementFunctor() {};

      void setGridIndex(size_t grid_index) override;
      size_t getNumGrids() override;

      /**
       * Evaluates the grids at necessary grid points and stores them for later
       * use. Needs to be called before a refinement step is applied.
       */
      void preComputeEvaluations() override;


    protected:
      std::vector<base::Grid*> grids;
      std::vector<base::DataVector*> alphas;

      size_t current_grid_index;
      size_t refinements_num;
      double threshold;
      bool level_penalize;

      bool pre_compute;
      std::vector<std::map<std::string, double>> pre_comp_evals;

    };

  }  // namespace base
}  // namespace sgpp

#endif /* GRIDPOINTBASEDREFINEMENTFUNCTOR_HPP */
