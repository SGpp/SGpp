// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertLinearToPrewavelet.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertPrewaveletToLinear.hpp>
#include <sgpp/base/operation/hash/OperationConvertPrewavelet.hpp>


#include <sgpp/base/algorithm/sweep.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationConvertPrewavelet::doConvertToLinear(
      DataVector& alpha) {

      ConvertPrewaveletToLinear func(this->storage);
      sweep<ConvertPrewaveletToLinear> s(func, this->storage);


      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }

    }

    void OperationConvertPrewavelet::doConvertFromLinear(DataVector& alpha) {
      ConvertLinearToPrewavelet func(this->storage, this->shadowstorage);
      sweep<ConvertLinearToPrewavelet> s(func, this->storage);

      for (size_t i = 0; i < this->storage->dim(); i++) {
        s.sweep1D(alpha, alpha, i);
      }
    }

  }
}