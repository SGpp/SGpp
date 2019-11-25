

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace datadriven {

struct VpHeapItem {
  VpHeapItem(size_t index, double distance) :
  index(index), distance(distance) {}
  size_t index;
  double distance;

  bool operator<(const VpHeapItem& o) const {
    return distance < o.distance;
  }
};
}  // namespace datadriven
}  // namespace sgpp
