// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPERATION_OPERATIONOPFACTORY_HPP
#define SGPP_OPTIMIZATION_OPERATION_OPERATIONOPFACTORY_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp>

namespace SGPP {
namespace op_factory {

/**
 * Creates a OperationMultipleHierarchisation for the given
 * SGPP::optimization grid.
 * Don't forget to delete the object after use.
 *
 * @param grid  sparse grid
 * @return      pointer to a OperationMultipleHierarchisation object
 *              for the grid
 */
optimization::OperationMultipleHierarchisation* createOperationMultipleHierarchisation(
    base::Grid& grid);
}  // namespace op_factory
}  // namespace SGPP

#endif /* SGPP_OPTIMIZATION_OPERATION_OPERATIONOPFACTORY_HPP */
