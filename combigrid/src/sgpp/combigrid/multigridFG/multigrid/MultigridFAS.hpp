// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIGRIDFAS_HPP_
#define MULTIGRIDFAS_HPP_

#include <sgpp/combigrid/multigridFG/interface/OperatorFG.hpp>

namespace combigrid {

  /** multigrid solver for a full grid. <br>
   * Actually it contains also Krylov space methods, and a variety of other methods
   * (hierarchical and non-hierarchical)
   * <br>
   * It uses an operator which contains the problem specific operators */
  class MultigridFAS {
    public:

      /** Ctor
       * @param op [IN] the operator which defines the problem
       * @param fg [IN] the full grid on which the problem should be solved
       * @param createHierarchy [IN] if a hierarchy of grids should be created */
      MultigridFAS(OperatorFG* op ,
                   FullGridD* fg ,
                   bool createHierarchy = true );

      /** */
      virtual ~MultigridFAS();

      /** solves the problem using the full approximation scheme
      * @param unknowns [IN/OUT] the unknown vector, which will be the initial solution
      * and later the output for the solution
      * @param errorTol [IN] the error tolerance for the solver */
      void solveFAS( std::vector<double>& unknowns , double errorTol);

    private:

      /** the hierarchy of grids for the multigrid,
       * the first one is the input grid, the rest will be created*/
      std::vector<FullGridD*> fullgrids_;

      /** operators for each full grid */
      std::vector<OperatorFG*> operators_;

      std::vector< std::vector<double>* > unknowns_;

      std::vector< std::vector<double>* > correction_;

      std::vector< std::vector<double>* > rhs_;

      std::vector< std::vector<double>* > rhs_tmp_;

      std::vector< std::vector<double>* > u_hH_;

      std::vector< std::vector<double>* > lh_;

      /** */
      int depth_;

      int nrGSPre_;

      int nrGSPost_;
  };

}

#endif /* MULTIGRIDFAS_HPP_ */