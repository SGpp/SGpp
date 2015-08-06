// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGULARIZATIONCONFIGURATION_HPP_
#define REGULARIZATIONCONFIGURATION_HPP_

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace pde {

    enum RegularizationType {
      Identity,
      Laplace
    };

    struct RegularizationConfiguration {
      RegularizationType regType_;
    };

  } // pde
} // SGPP

#endif /* REGULARIZATIONCONFIGURATION_HPP_ */
