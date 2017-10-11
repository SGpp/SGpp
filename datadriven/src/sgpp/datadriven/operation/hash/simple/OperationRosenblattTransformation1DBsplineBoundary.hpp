// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONROSENBLATTTRANSFORMATION1DBSPLINEBOUNDARY_HPP
#define OPERATIONROSENBLATTTRANSFORMATION1DBSPLINEBOUNDARY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
class OperationRosenblattTransformation1DBsplineBoundary : public OperationTransformation1D {
 protected:
  base::Grid* grid;

 public:
  explicit OperationRosenblattTransformation1DBsplineBoundary(base::Grid* grid);
  virtual ~OperationRosenblattTransformation1DBsplineBoundary();

  /**
   * Rosenblatt Transformation 1D
   * @param alpha1d
   * @param coord1d
   * @return
   */
  double doTransformation1D(base::DataVector* alpha1d, double coord1d);
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONROSENBLATTTRANSFORMATION1DBSPLINEBOUNDARY_HPP */
