// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/json/JSON.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

class DataMiningConfiguration {
 public:
  DataMiningConfiguration();
};

}  // namespace datadriven
}  // namespace sgpp
