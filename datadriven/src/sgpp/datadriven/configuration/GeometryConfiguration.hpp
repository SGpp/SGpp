/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GeometryConfiguration.hpp
 *
 *  Created on: Jan 15, 2019
 *      Author: jamal
 */

#pragma once

#include <sgpp/globaldef.hpp>
#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {

enum class StencilType { DN, None };

/*
 * Struct that stores information to geometry aware sparse grids
 */
struct GeometryConfiguration{
  /*
   * Stencil for geometric relation
   */
  StencilType stencilType;

  /*
   * resolution of image/video e.g 28x28
   */
  std::vector<int64_t> dim;
};


}  // namespace datadriven
}  // namespace sgpp
