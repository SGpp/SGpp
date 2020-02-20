// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/LevelIndexTypes.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

using sgpp::base::index_t;
using sgpp::base::level_t;

/// level multi-index
typedef std::vector<level_t> LevelVector;
/// index multi-index
typedef std::vector<index_t> IndexVector;


/**
 * @brief get a hypercube of all level vectors between and including minLevel and maxLevel
 *
 * @param maxLevel                    the maximum level vector
 * @param minLevel                    the minimum level vector, must not be larger than maxLevel in
 *                                      any dimension
 * @return std::vector<LevelVector>   a lexicographically ordered vector of level vectors in the
 *                                      hypercube //TODO is this true?
 */
std::vector<LevelVector> hyperCubeOfLevelVectors(const LevelVector& maxLevel,
                                                 const LevelVector& minLevel);

}  // namespace combigrid
}  // namespace sgpp
