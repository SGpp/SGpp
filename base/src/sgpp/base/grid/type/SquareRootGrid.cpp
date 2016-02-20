// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>

#include <sgpp/base/grid/generation/SquareRootGridGenerator.hpp>
// #include <sgpp/base/operation/hash/OperationEvalLinearBoundary.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

SquareRootGrid::SquareRootGrid(std::istream& istr) :
  Grid(istr) {
}

SquareRootGrid::SquareRootGrid(size_t dim) :
  Grid(dim) {
}

SquareRootGrid::SquareRootGrid(BoundingBox& BB) :
  Grid(BB) {
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
std::unique_ptr<GridGenerator> SquareRootGrid::createGridGenerator() {
  return std::unique_ptr<GridGenerator>(new SquareRootGridGenerator(*this->storage));
}
// OperationHierarchisation* SquareRootGrid::createOperationHierarchisation()
// {
//   return new OperationHierarchisationLinearBoundary(*this->storage);
// }
// OperationEval* SquareRootGrid::createOperationEval()
// {
//   return new OperationEvalLinearBoundary(*this->storage);
// }

// OperationConvert* SquareRootGrid::createOperationConvert()
// {
//   throw factory_exception("Unsupported operation");
// }

}  // namespace base
}  // namespace SGPP
