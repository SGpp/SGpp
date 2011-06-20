/*
 * CombiArbitraryScheme.hpp
 *
 *  Created on: Apr 28, 2011
 *      Author: kowitz_local
 */

#ifndef COMBIARBITRARYSCHEME_HPP_
#define COMBIARBITRARYSCHEME_HPP_

#include "CombiSchemeBasis.hpp"
#include "combigrid/utils/CombigridLevelVector.hpp"

namespace combigrid {
/** Combischeme created with from an arbitrary active set of full grids. Missing
 * subgrids for the creation of a valid combi solution are automatically added.
 */
class CombiArbitraryScheme: public combigrid::CombiSchemeBasis {
public:
	/** Initializes the arbitrary scheme from a set of active set level vectors*/
	CombiArbitraryScheme(std::vector<std::vector<int> > level_vectors);
	/** Initializes the arbitrary scheme from a CombigridLevelVector as acitve set.*/
	CombiArbitraryScheme(CombigridLevelVector in);

};

}

#endif /* COMBIARBITRARYSCHEME_HPP_ */
