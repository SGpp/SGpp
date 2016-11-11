// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <sgpp/combigrid/MultiFunction.hpp>

namespace sgpp {
namespace combigrid {

/**
 * This is a helper class used internally to wrap python functions into MultiFunction objects.
 */
class MultiFunctionDirector {
 public:
  virtual ~MultiFunctionDirector();

  virtual double eval(base::DataVector const &vec) = 0;

  virtual MultiFunction toMultiFunction();
};

} /* namespace combigrid */
} /* namespace sgpp */
