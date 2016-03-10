// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTRANSFORMATION1D_HPP
#define OPERATIONTRANSFORMATION1D_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>
#include <cstring>

namespace sgpp {
namespace datadriven {

/**
 * Sample 1D Probability Density Function
 */

class OperationTransformation1D {
 public:
  OperationTransformation1D() {}
  virtual ~OperationTransformation1D() {}

  /**
   * Transform 1d
   * @param alpha1d
   * @param coord1d
   * @return
   */
  virtual double doTransformation1D(base::DataVector* alpha1d, double coord1d) = 0;
};

}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONTRANSFORMATION1D_HPP */
