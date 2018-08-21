// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/FundamentalNotAKnotSplineBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

FundamentalNotAKnotSplineBoundaryGrid::FundamentalNotAKnotSplineBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SFundamentalNotAKnotSplineBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

FundamentalNotAKnotSplineBoundaryGrid::FundamentalNotAKnotSplineBoundaryGrid(
    size_t dim, size_t degree, level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SFundamentalNotAKnotSplineBase(degree)),
      boundaryLevel(boundaryLevel) {}

FundamentalNotAKnotSplineBoundaryGrid::~FundamentalNotAKnotSplineBoundaryGrid() {}

sgpp::base::GridType FundamentalNotAKnotSplineBoundaryGrid::getType() {
  return sgpp::base::GridType::FundamentalNotAKnotSplineBoundary;
}

SBasis& FundamentalNotAKnotSplineBoundaryGrid::getBasis() { return *basis_; }

size_t FundamentalNotAKnotSplineBoundaryGrid::getDegree() { return this->degree; }

Grid* FundamentalNotAKnotSplineBoundaryGrid::unserialize(std::istream& istr) {
  return new FundamentalNotAKnotSplineBoundaryGrid(istr);
}

void FundamentalNotAKnotSplineBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& FundamentalNotAKnotSplineBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
