/*
 * CombiArbitraryScheme.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: kowitz_local
 */

#include "CombiArbitraryScheme.hpp"

using namespace combigrid;

CombiArbitraryScheme::CombiArbitraryScheme(
		std::vector<std::vector<int> > level_vectors) :
	CombiSchemeBasis(level_vectors[0].size(), level_vectors[0]) {
	CombigridLevelVector buffer = CombigridLevelVector::getCombiLevels(
			level_vectors);
	levels_vector_ = buffer.getLevelVec();
	cofficients_ = buffer.getCoef();
}

