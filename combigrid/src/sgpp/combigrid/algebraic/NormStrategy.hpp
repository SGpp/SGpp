// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

namespace sgpp {
namespace combigrid {

template <typename V>
class NormStrategy {
 public:
  NormStrategy(){};
  virtual ~NormStrategy(){};

  double norm(V vector) { return vector.norm(); }
};

} /* namespace combigrid */
} /* namespace sgpp */
