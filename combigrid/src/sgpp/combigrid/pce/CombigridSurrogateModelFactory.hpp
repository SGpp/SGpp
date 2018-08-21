// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

namespace sgpp {
namespace combigrid {

/**
 *
 * @param config
 * @return
 */
std::shared_ptr<CombigridSurrogateModel> createCombigridSurrogateModel(
    CombigridSurrogateModelConfiguration& config);

} /* namespace combigrid */
} /* namespace sgpp */
