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

class CombiArbitraryScheme: public combigrid::CombiSchemeBasis {
public:
	CombiArbitraryScheme(std::vector<std::vector<int> > level_vectors);
	CombiArbitraryScheme(CombigridLevelVector in);

};

}

#endif /* COMBIARBITRARYSCHEME_HPP_ */
