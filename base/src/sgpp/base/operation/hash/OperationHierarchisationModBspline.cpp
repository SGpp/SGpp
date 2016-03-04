// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationModBspline.hpp>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void OperationHierarchisationModBspline::doHierarchisation(
  DataVector& node_values) {
  throw operation_exception(
    "This operation is not implemented, yet! Sorry ;-)");
}

void OperationHierarchisationModBspline::doDehierarchisation(
  DataVector& alpha) {
  throw operation_exception(
    "This operation is not implemented, yet! Sorry ;-)");
}

}  // namespace base
}  // namespace sgpp
