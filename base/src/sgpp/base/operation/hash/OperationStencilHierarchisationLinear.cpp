// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationStencilHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/StencilHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/StencilDehierarchisationLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

void OperationStencilHierarchisationLinear::doHierarchisation(
  DataVector& node_values) {
  surplusStencil.clear();
  neighborStencil.clear();
  weightStencil.clear();
  StencilHierarchisationLinear func(this->storage, surplusStencil,
                                    neighborStencil, weightStencil);
  sweep<StencilHierarchisationLinear> s(func, this->storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage->dim(); i++) {
    s.sweep1D(node_values, node_values, i);
  }
}

void OperationStencilHierarchisationLinear::doDehierarchisation(
  DataVector& alpha) {
  surplusStencil.clear();
  neighborStencil.clear();
  weightStencil.clear();
  StencilDehierarchisationLinear func(this->storage, surplusStencil,
                                      neighborStencil, weightStencil);
  sweep<StencilDehierarchisationLinear> s(func, this->storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage->dim(); i++) {
    s.sweep1D(alpha, alpha, i);
  }
}

}  // namespace base
}  // namespace SGPP
