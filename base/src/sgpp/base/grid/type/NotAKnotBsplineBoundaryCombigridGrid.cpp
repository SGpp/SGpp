// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/NotAKnotBsplineBoundaryCombigridGrid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NotAKnotBsplineBoundaryCombigridGrid::NotAKnotBsplineBoundaryCombigridGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SNotAKnotBsplineBoundaryCombigridBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

NotAKnotBsplineBoundaryCombigridGrid::NotAKnotBsplineBoundaryCombigridGrid(size_t dim,
                                                                           size_t degree,
                                                                           level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SNotAKnotBsplineBoundaryCombigridBase(degree)),
      boundaryLevel(boundaryLevel) {}

NotAKnotBsplineBoundaryCombigridGrid::~NotAKnotBsplineBoundaryCombigridGrid() {}

sgpp::base::GridType NotAKnotBsplineBoundaryCombigridGrid::getType() {
  return sgpp::base::GridType::NotAKnotBsplineBoundary;
}

const SBasis& NotAKnotBsplineBoundaryCombigridGrid::getBasis() { return *basis_; }

size_t NotAKnotBsplineBoundaryCombigridGrid::getDegree() { return this->degree; }

Grid* NotAKnotBsplineBoundaryCombigridGrid::unserialize(std::istream& istr) {
  return new NotAKnotBsplineBoundaryCombigridGrid(istr);
}

void NotAKnotBsplineBoundaryCombigridGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NotAKnotBsplineBoundaryCombigridGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
