// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModPolyClenshawCurtisGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

ModPolyClenshawCurtisGrid::ModPolyClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SPolyModifiedClenshawCurtisBase(degree));
}

ModPolyClenshawCurtisGrid::ModPolyClenshawCurtisGrid(size_t dim, size_t degree)
    : Grid(dim),
      generator(storage),
      degree(degree),
      basis_(new SPolyModifiedClenshawCurtisBase(degree)) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

ModPolyClenshawCurtisGrid::~ModPolyClenshawCurtisGrid() {}

sgpp::base::GridType ModPolyClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::ModPolyClenshawCurtis;
}

SBasis& ModPolyClenshawCurtisGrid::getBasis() { return *basis_; }

size_t ModPolyClenshawCurtisGrid::getDegree() { return this->degree; }

Grid* ModPolyClenshawCurtisGrid::unserialize(std::istream& istr) {
  return new ModPolyClenshawCurtisGrid(istr);
}

void ModPolyClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModPolyClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
