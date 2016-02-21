// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>


namespace SGPP {
namespace base {

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream& istr) :
  Grid(istr),
  boundaryLevel(0) {
  istr >> boundaryLevel;
}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim,
    level_t boundaryLevel) :
  Grid(dim),
  boundaryLevel(boundaryLevel) {
}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(BoundingBox& BB,
    level_t boundaryLevel) :
  Grid(BB),
  boundaryLevel(boundaryLevel) {
}

LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid() {
}

SGPP::base::GridType LinearClenshawCurtisGrid::getType() {
  return SGPP::base::GridType::LinearClenshawCurtis;
}

const SBasis& LinearClenshawCurtisGrid::getBasis() {
  static SLinearClenshawCurtisBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearClenshawCurtisGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearClenshawCurtisGrid(istr));
}

void LinearClenshawCurtisGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
std::unique_ptr<GridGenerator> LinearClenshawCurtisGrid::createGridGenerator() {
  return std::unique_ptr<GridGenerator>(new BoundaryGridGenerator(storage, boundaryLevel));
}


}  // namespace base
}  // namespace SGPP
