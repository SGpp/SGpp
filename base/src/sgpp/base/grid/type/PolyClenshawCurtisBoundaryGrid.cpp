// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/PolyClenshawCurtisBoundaryGrid.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <vector>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp>

namespace sgpp {
namespace base {

PolyClenshawCurtisBoundaryGrid::PolyClenshawCurtisBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), boundaryLevel(0), degree(0) {
  istr >> degree;
  istr >> boundaryLevel;

  generator.setBoundaryLevel(boundaryLevel);
  basis_.reset(new SPolyClenshawCurtisBoundaryBase(degree));
}

PolyClenshawCurtisBoundaryGrid::PolyClenshawCurtisBoundaryGrid(size_t dim, size_t degree,
                                                               level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      boundaryLevel(boundaryLevel),
      degree(degree),
      basis_(new SPolyClenshawCurtisBoundaryBase(degree)) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

PolyClenshawCurtisBoundaryGrid::~PolyClenshawCurtisBoundaryGrid() {}

size_t PolyClenshawCurtisBoundaryGrid::getDegree() const { return this->degree; }

sgpp::base::GridType PolyClenshawCurtisBoundaryGrid::getType() {
  return sgpp::base::GridType::PolyClenshawCurtisBoundary;
}

SBasis& PolyClenshawCurtisBoundaryGrid::getBasis() { return *basis_; }

Grid* PolyClenshawCurtisBoundaryGrid::unserialize(std::istream& istr) {
  return new PolyClenshawCurtisBoundaryGrid(istr);
}

void PolyClenshawCurtisBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& PolyClenshawCurtisBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
