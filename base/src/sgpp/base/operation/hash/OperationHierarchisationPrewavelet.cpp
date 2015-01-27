// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationPrewavelet.hpp>
#include <sgpp/base/basis/prewavelet/algorithm_sweep/ConvertLinearToPrewavelet.hpp>
#include <sgpp/base/basis/prewavelet/algorithm_sweep/ConvertPrewaveletToLinear.hpp>
#include <sgpp/base/basis/linear/noboundary/algorithm_sweep/HierarchisationLinear.hpp>
#include <sgpp/base/basis/linear/noboundary/algorithm_sweep/DehierarchisationLinear.hpp>


#include <sgpp/base/algorithm/sweep.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationHierarchisationPrewavelet::doHierarchisation(
      DataVector& node_values) {
      /*
       * Hierarchisation on prewavelets require a hierarchisation on a normal
       * linear grid, afterwards they are converted into a prewavelet basis
       */

      HierarchisationLinear func(this->storage);
      sweep<HierarchisationLinear> s(func, this->storage);

      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(node_values, node_values, i);
      }

      ConvertLinearToPrewavelet func2(this->storage, this->shadowStorage);
      sweep<ConvertLinearToPrewavelet> s2(func2, this->storage);

      for (size_t i = 0; i < this->storage->dim(); i++) {
        s2.sweep1D(node_values, node_values, i);
      }

    }

    void OperationHierarchisationPrewavelet::doDehierarchisation(DataVector& alpha) {
      ConvertPrewaveletToLinear func(this->storage);
      sweep<ConvertPrewaveletToLinear> s(func, this->storage);

      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }

      DehierarchisationLinear func2(this->storage);
      sweep<DehierarchisationLinear> s2(func2, this->storage);

      for (size_t i = 0; i < this->storage->dim(); i++) {
        s2.sweep1D(alpha, alpha, (this->storage->dim() - (i + 1)));
      }

    }

    void OperationHierarchisationPrewavelet::expandGrid() {
      for (size_t i = 0; i < shadowStorage->size(); i++) {
        (*shadowStorage->get(i)).toString(std::cout);
        this->storage->insert(*shadowStorage->get(i));

        if ((*shadowStorage->get(i)).isLeaf())
          std::cout << "is Leaf : " << std::endl;
        else
          std::cout << "nööö" << std::endl;
      }
    }

    void OperationHierarchisationPrewavelet::shrinkGrid() {
      for (size_t i = 0; i < shadowStorage->size(); i++) {
        this->storage->deleteLast();
      }
    }

  }
}