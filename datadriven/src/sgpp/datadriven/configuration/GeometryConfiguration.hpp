// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

enum class StencilType {
  DirectNeighbour,
  AllHierarchicalParent,
  NextHierarchicalParent,
  Block,
  None
};

struct StencilConfiguration {
  /*
   * Stenciltype to be used
   */
  StencilType stencilType_;
  /*
   * Determines on which layers the stencils should be applied
   */
  std::vector<size_t> applyOnLayers_;
  /*
   * Index of the color channel for this specific stencil
   * -1 if no color channels available
   */
  int64_t colorIndex_;
  /*
   * Blocklenght for Blockstencil
   */
  size_t blockLenght_;
};

/*
 * Struct that stores information to geometry aware sparse grids
 */
struct GeometryConfiguration {
  /*
   * Stencil for geometric relation
   */
  std::vector<StencilConfiguration> stencils_;

  /*
   * resolution of image/video e.g 28x28
   */
  std::vector<std::vector<int64_t>> dim_;
};

}  // namespace datadriven
}  // namespace sgpp
