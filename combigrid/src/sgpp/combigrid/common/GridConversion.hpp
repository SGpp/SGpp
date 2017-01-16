// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

namespace sgpp {
namespace combigrid {

std::shared_ptr<TreeStorage<uint8_t>> allLevelsInHashGridStorage(base::HashGridStorage *storage);
base::HashGridStorage *toHashGridStorage(std::shared_ptr<TreeStorage<uint8_t>> levelStructure);

} /* namespace combigrid */
} /* namespace sgpp */
