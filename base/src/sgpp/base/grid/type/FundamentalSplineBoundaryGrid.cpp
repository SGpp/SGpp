// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/FundamentalSplineBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

FundamentalSplineBoundaryGrid::FundamentalSplineBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SFundamentalSplineBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

FundamentalSplineBoundaryGrid::FundamentalSplineBoundaryGrid(
    size_t dim, size_t degree, level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SFundamentalSplineBase(degree)),
      boundaryLevel(boundaryLevel) {}

FundamentalSplineBoundaryGrid::~FundamentalSplineBoundaryGrid() {}

sgpp::base::GridType FundamentalSplineBoundaryGrid::getType() {
  return sgpp::base::GridType::FundamentalSplineBoundary;
}

SBasis& FundamentalSplineBoundaryGrid::getBasis() { return *basis_; }

size_t FundamentalSplineBoundaryGrid::getDegree() { return this->degree; }

Grid* FundamentalSplineBoundaryGrid::unserialize(std::istream& istr) {
  return new FundamentalSplineBoundaryGrid(istr);
}

void FundamentalSplineBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& FundamentalSplineBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
