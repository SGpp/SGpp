/*
 * Multigrid.hpp
 *
 *  Created on: May 16, 2011
 *      Author: benk
 */

#ifndef MULTIGRID_HPP_
#define MULTIGRID_HPP_

#include "solver/multigridFG/interface/OperatorFG.hpp"

namespace combigrid {

/** multigrid solver for a full grid. <br>
 * It uses an operator which contains the problem specific operators */
class Multigrid {
public:

	/** Ctor
	 * @param op [IN] the operator which defines the problem
	 * @param fg [IN] the full grid on which the problem should be solved
	 * @param createHierarchy [IN] if a hierarchy of grids should be created */
	Multigrid(OperatorFG* op ,
			  FullGridD* fg ,
			  bool createHierarchy = true );

	/** */
	virtual ~Multigrid();

	/** solves the problem using the correction scheme
	 * @param unknowns [IN/OUT] the unknown vector, which will be the initial solution
	 * and later the output for the solution
	 * @param errorTol [IN] the error tolerance for the solver*/
    void solveCS( std::vector<double>& unknowns , double errorTol);

    /** solves the problem using the full approximation scheme
	 * @param unknowns [IN/OUT] the unknown vector, which will be the initial solution
	 * and later the output for the solution
	 * @param errorTol [IN] the error tolerance for the solver */
    void solveFAS( std::vector<double>& unknowns , double errorTol);

	/** solves the problem using just the smoothing
	 * @param unknowns [IN/OUT] the unknown vector, which will be the initial solution
	 * and later the output for the solution
	 * @param errorTol [IN] the error tolerance for the solver*/
    void solveSmoothing( std::vector<double>& unknowns , double errorTol);

private:

    /** the hierarchy of grids for the multigrid,
     * the first one is the input grid, the rest will be created*/
    std::vector<FullGridD*> fullgrids_;

    /** operators for each full grid */
    std::vector<OperatorFG*> operators_;

    std::vector< std::vector<double>* > unknowns_;

    std::vector< std::vector<double>* > correction_;

    std::vector< std::vector<double>* > rhs_;

    /** */
    int depth_;

    int nrGSPre_;

    int nrGSPost_;
};

}

#endif /* MULTIGRID_HPP_ */
