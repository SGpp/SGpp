// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

PolyClenshawCurtisGrid::PolyClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SPolyClenshawCurtisBase(degree));
}

PolyClenshawCurtisGrid::PolyClenshawCurtisGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SPolyClenshawCurtisBase(degree)) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

PolyClenshawCurtisGrid::~PolyClenshawCurtisGrid() {}

sgpp::base::GridType PolyClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::PolyClenshawCurtis;
}

SBasis& PolyClenshawCurtisGrid::getBasis() { return *basis_; }

size_t PolyClenshawCurtisGrid::getDegree() { return this->degree; }

Grid* PolyClenshawCurtisGrid::unserialize(std::istream& istr) {
  return new PolyClenshawCurtisGrid(istr);
}

void PolyClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& PolyClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
