// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

using sgpp::base::level_t;
using sgpp::base::index_t;
typedef std::vector<level_t> LevelVector;
typedef std::vector<index_t> IndexVector;

}  // namespace combigrid
}  // namespace sgpp
