/*
 * RunTikhonov.hpp
 *
 *  Created on: May 30, 2011
 *      Author: benk
 */

#ifndef RUNTIKHONOV_HPP_
#define RUNTIKHONOV_HPP_

#include "combigrid.hpp"
#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/combigrid/AbstractCombiGrid.hpp"
#include "combigrid/domain/CombiGridDomain.hpp"

namespace combigrid {

/** class to run the test cases for tikhonov regularization <br>
 * The input */
class RunTikhonov {
public:

	/** empty Ctor*/
	RunTikhonov() {;}

	/** empty Dtor*/
	virtual ~RunTikhonov() {;}

	/** static function to run the Tikhonov reguralization on a full grid
	 * @param domain [IN] the domain for the full grid on which the regularization will be done
	 * @param levels [IN] level vector
	 * @param lambda [IN] the lambda factor
	 * @param Xfile  [IN] file path with the X coordinates
	 * @param Yfile  [IN] file with the Y values
	 * @return a full grid with the solution */
	static FullGridD* computeFGTikhonov(
			GridDomain& domain,
			const std::vector<int>& levels,
			double lambda,
			const std::string& Xfile ,
			const std::string& YFile );

private :

	/** read in the input values for the regularization*/
	static void readInInput(const std::string& XfileS , const std::string& YFileS ,
			         int& dimensions ,  int& nrPoints ,
			         std::vector<double>& XCoords, std::vector<double>& YPoint);
};

}

#endif /* RUNTIKHONOV_HPP_ */
