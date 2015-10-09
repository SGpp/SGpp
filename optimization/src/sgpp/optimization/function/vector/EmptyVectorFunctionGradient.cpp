// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunctionGradient.hpp>

namespace SGPP {
  namespace optimization {

    WrapperVectorFunctionGradient emptyVectorFunctionGradient(
      0, 0, [](const base::DataVector& x,
               base::DataVector& value,
    base::DataMatrix& gradient) {});

  }
}
