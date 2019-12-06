// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/operation/OperationPoleNodalisationLinear.hpp>

namespace sgpp {
namespace combigrid {

OperationPoleNodalisationLinear::OperationPoleNodalisationLinear() {
}

OperationPoleNodalisationLinear::~OperationPoleNodalisationLinear() {
}

void OperationPoleNodalisationLinear::apply(base::DataVector& values, size_t start, size_t step,
    size_t count, level_t level, bool hasBoundary) {
  // do nothing, as nodal coefficients equal values
}

}  // namespace combigrid
}  // namespace sgpp
