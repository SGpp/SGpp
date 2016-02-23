// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

SquareRootGrid::SquareRootGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

SquareRootGrid::SquareRootGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

SquareRootGrid::SquareRootGrid(BoundingBox& BB) :
  Grid(BB),
  generator(storage) {
}

SquareRootGrid::~SquareRootGrid() {
}

SGPP::base::GridType SquareRootGrid::getType() {
  return SGPP::base::GridType::SquareRoot;
}

const SBasis& SquareRootGrid::getBasis() {
  static SLinearBoundaryBase basis;
  return basis;
}

std::unique_ptr<Grid> SquareRootGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new SquareRootGrid(istr));
}
/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& SquareRootGrid::getGenerator() {
  return generator;
}
// OperationHierarchisation* SquareRootGrid::createOperationHierarchisation()
// {
//   return new OperationHierarchisationLinearBoundary(storage);
// }
// OperationEval* SquareRootGrid::createOperationEval()
// {
//   return new OperationEvalLinearBoundary(storage);
// }

// OperationConvert* SquareRootGrid::createOperationConvert()
// {
//   throw factory_exception("Unsupported operation");
// }

}  // namespace base
}  // namespace SGPP
