// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "CombiArbitraryScheme.hpp"

combigrid::CombiArbitraryScheme::CombiArbitraryScheme(std::vector < std::vector <
    int > > level_vectors) :
  combigrid::CombiSchemeBasis(static_cast<int>(level_vectors[0].size()), level_vectors[0]) {
  combigrid::CombigridLevelVector buffer =
    combigrid::CombigridLevelVector::getCombiLevels(level_vectors);
  levels_vector_ = buffer.getLevelVec();
  cofficients_ = buffer.getCoef();
}

combigrid::CombiArbitraryScheme::CombiArbitraryScheme(
  combigrid::CombigridLevelVector in) :
  combigrid::CombiSchemeBasis(in.getDim(), in.getLevelVec()[0]) {
  levels_vector_ = in.getLevelVec();
  cofficients_ = in.getCoef();
}