// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/FundamentalNakSplineBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

FundamentalNakSplineBoundaryGrid::FundamentalNakSplineBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SFundamentalNakSplineBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

FundamentalNakSplineBoundaryGrid::FundamentalNakSplineBoundaryGrid(
    size_t dim, size_t degree, level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SFundamentalNakSplineBase(degree)),
      boundaryLevel(boundaryLevel) {}

FundamentalNakSplineBoundaryGrid::~FundamentalNakSplineBoundaryGrid() {}

sgpp::base::GridType FundamentalNakSplineBoundaryGrid::getType() {
  return sgpp::base::GridType::FundamentalNakSplineBoundary;
}

SBasis& FundamentalNakSplineBoundaryGrid::getBasis() { return *basis_; }

size_t FundamentalNakSplineBoundaryGrid::getDegree() { return this->degree; }

Grid* FundamentalNakSplineBoundaryGrid::unserialize(std::istream& istr) {
  return new FundamentalNakSplineBoundaryGrid(istr);
}

void FundamentalNakSplineBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& FundamentalNakSplineBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
