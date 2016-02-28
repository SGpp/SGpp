// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/json/JSON.hpp>

#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <vector>
#include <map>
#include <string>

namespace SGPP {
namespace datadriven {

class DataMiningConfiguration : public json::JSON {
 public:
  DataMiningConfiguration();

  explicit DataMiningConfiguration(const std::string& fileName);

  virtual DataMiningConfiguration* clone();

  base::GridType stringToGridType(std::string& gridType);
  RegularizationType stringToRegularizationType(std::string& regularizationType);
  solver::SLESolverType stringToSolverType(std::string& solverType);
};

}  // namespace datadriven
}  // namespace SGPP
