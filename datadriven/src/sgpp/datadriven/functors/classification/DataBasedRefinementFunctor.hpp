// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#ifndef DATABASEDREFINEMNTFUNCTOR_HPP
#define DATABASEDREFINEMNTFUNCTOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>

#include <vector>


namespace sgpp {
namespace datadriven {

  /**
   * Data based refinement uses data points to find refinement candidates.
   * For the given data sets class-intersection sets H_k are computed. A
   * grid points is included in H_k if for at least on class l
   * PDF_k(point) > coeff_a_k * mu_k AND PDF_l(point) > coeff_a_l * mu_l.
   * To determine the score of a grid point, the number of data points
   * from H_k within the support of this grid point is taken
   */
  class DataBasedRefinementFunctor : public MultiGridRefinementFunctor {
    public:
      /**
       * Constructor.
       *
       * @param grids Vector of grids. current_grid_index specifies the grid to be refined
       * @param alphas Vector of surpluses related to the grids
       * @param data The data used to compute the sets H_k
       * @param targets The classes for data
       * @param refinements_num Maximum number of refinements done
       * @param level_penalize If a level penalizing is multiplied to the score (2^{|l|_1})
       * @param coeff_a Scaling coefficients for the computation of H_k. Per default 1.0
       * @param threshold Threshold for refinement scores
       */
      DataBasedRefinementFunctor(std::vector<base::Grid*> grids,
                                 std::vector<base::DataVector*> alphas,
                                 base::DataMatrix* data,
                                 base::DataVector* targets,
                                 size_t refinements_num = 1,
                                 bool level_penalize = false,
                                 std::vector<double> coeff_a =
                                 std::vector<double>(),
                                 double threshold = 0.0);

      double operator()(base::GridStorage& storage,
                        size_t seq) const override;
      double start() const override;
      size_t getRefinementsNum() const override;
      double getRefinementThreshold() const override;
      virtual ~DataBasedRefinementFunctor() {}

      void setGridIndex(size_t grid_index) override;
      size_t getNumGrids() override;

      // Sets the data/target pointer
      void setData(base::DataMatrix* data, base::DataVector* targets);

      // (re-)computes all H_k based on data and targets
      void computeH();

      // Get the set Hk (for plotting/debugging)
      base::DataMatrix& getHk(size_t index);


    protected:
      std::vector<base::Grid*> grids;
      std::vector<base::DataVector*> alphas;
      base::DataMatrix evals;

      base::DataMatrix* data;
      base::DataVector* targets;

      // The set H containing all H_k
      std::vector<base::DataMatrix> h;

      // Means of the PDF for all classes
      std::vector<double> means;

      // Scaling coefficient of the means
      std::vector<double> coeff_a;

      size_t current_grid_index;
      size_t refinements_num;
      double threshold;
      bool level_penalize;


      // cl_indx is a index into classes rather than the label itself
      void computeHkl(base::DataMatrix& inters,
                      size_t cl_ind1,
                      size_t cl_ind2);
      bool isWithinSupport(base::HashGridPoint& gp,
                           base::DataVector& point) const;
    };

}  // namespace base
}  // namespace sgpp

#endif /* DATABASEDREFINEMENTFUNCTOR_HPP */
