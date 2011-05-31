/*
 * RunPoisson.hpp
 *
 *  Created on: May 30, 2011
 *      Author: benk
 */

#ifndef RunPoisson_HPP_
#define RunPoisson_HPP_

#include "combigrid.hpp"
#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/combigrid/AbstractCombiGrid.hpp"
#include "combigrid/domain/CombiGridDomain.hpp"
#include "solver/multigridFG/operators/PoissonOperator.hpp"

namespace combigrid {

/** class to run the test cases for tikhonov regularization <br>
 * The input */
class RunPoisson {
public:

	/** empty Ctor*/
	RunPoisson() {;}

	/** empty Dtor*/
	virtual ~RunPoisson() {;}

	/** static function to run the Poisson problem on a full grid
	 * @param domain [IN] the domain for the full grid
	 * @param levels [IN] level vector
	 * @param sigma [IN] diffusion coefficients in each direction
	 * @param startValue [IN] the value of all unknowns at the beginning
	 * @param callbackRHS  [IN] right hand side call back
	 * */
	static FullGridD* computeFGPoisson(
			GridDomain& domain,
			const std::vector<int>& levels,
			const std::vector<double>& sigma ,
			double startValue ,
			const CallBackRHS* callbackRHS );

private :

};

}

#endif /* RunPoisson_HPP_ */
