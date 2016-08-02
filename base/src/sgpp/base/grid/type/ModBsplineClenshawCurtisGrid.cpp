// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

ModBsplineClenshawCurtisGrid::ModBsplineClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SBsplineModifiedClenshawCurtisBase(degree));
}

ModBsplineClenshawCurtisGrid::ModBsplineClenshawCurtisGrid(size_t dim, size_t degree)
    : Grid(dim),
      generator(storage),
      degree(degree),
      basis_(new SBsplineModifiedClenshawCurtisBase(degree)) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

ModBsplineClenshawCurtisGrid::~ModBsplineClenshawCurtisGrid() {}

sgpp::base::GridType ModBsplineClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::ModBsplineClenshawCurtis;
}

const SBasis& ModBsplineClenshawCurtisGrid::getBasis() { return *basis_; }

size_t ModBsplineClenshawCurtisGrid::getDegree() { return this->degree; }

std::unique_ptr<Grid> ModBsplineClenshawCurtisGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new ModBsplineClenshawCurtisGrid(istr));
}

void ModBsplineClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModBsplineClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
