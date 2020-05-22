// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Configuration structure for the offline datamatrix decomposition database
 */
struct DatabaseConfiguration {
  /**
   * Filepath to the database
   */
  std::string filePath_ = "";
};

}  // namespace datadriven
}  // namespace sgpp
