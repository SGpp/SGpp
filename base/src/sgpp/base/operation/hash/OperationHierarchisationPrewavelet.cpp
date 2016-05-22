// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationPrewavelet.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertLinearToPrewavelet.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertPrewaveletToLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationLinear.hpp>


#include <sgpp/base/algorithm/sweep.hpp>



#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void OperationHierarchisationPrewavelet::doHierarchisation(
  DataVector& node_values) {
  /*
   * Hierarchisation on prewavelets require a hierarchisation on a normal
   * linear grid, afterwards they are converted into a prewavelet basis
   */

  HierarchisationLinear func(storage);
  sweep<HierarchisationLinear> s(func, storage);

  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(node_values, node_values, i);
  }

  ConvertLinearToPrewavelet func2(storage, shadowStorage);
  sweep<ConvertLinearToPrewavelet> s2(func2, storage);

  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s2.sweep1D(node_values, node_values, i);
  }
}

void OperationHierarchisationPrewavelet::doDehierarchisation(
  DataVector& alpha) {
  ConvertPrewaveletToLinear func(storage);
  sweep<ConvertPrewaveletToLinear> s(func, storage);

  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D(alpha, alpha, i);
  }

  DehierarchisationLinear func2(storage);
  sweep<DehierarchisationLinear> s2(func2, storage);

  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s2.sweep1D(alpha, alpha, (this->storage.getDimension() - (i + 1)));
  }
}

void OperationHierarchisationPrewavelet::expandGrid() {
  for (size_t i = 0; i < shadowStorage.getSize(); i++) {
    shadowStorage.getGridPoint(i).toString(std::cout);
    this->storage.insert(shadowStorage.getGridPoint(i));

    if (shadowStorage.getGridPoint(i).isLeaf())
      std::cout << "is Leaf : " << std::endl;
    else
      std::cout << "nooo" << std::endl;
  }
}

void OperationHierarchisationPrewavelet::shrinkGrid() {
  for (size_t i = 0; i < shadowStorage.getSize(); i++) {
    this->storage.deleteLast();
  }
}

}  // namespace base
}  // namespace sgpp
