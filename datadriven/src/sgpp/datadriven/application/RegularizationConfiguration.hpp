// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGULARIZATIONCONFIGURATION_HPP_
#define REGULARIZATIONCONFIGURATION_HPP_

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace datadriven {

enum class RegularizationType { Identity, Laplace, Diagonal, Lasso, ElasticNet };

struct RegularizationConfiguration {
  RegularizationType regType_;
  double lambda_;
  double l1Ratio_;
  double exponentBase_;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* REGULARIZATIONCONFIGURATION_HPP_ */
