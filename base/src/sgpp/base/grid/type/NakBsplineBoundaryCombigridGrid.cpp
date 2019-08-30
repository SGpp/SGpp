// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/NakBsplineBoundaryCombigridGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NakBsplineBoundaryCombigridGrid::NakBsplineBoundaryCombigridGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SNakBsplineBoundaryCombigridBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

NakBsplineBoundaryCombigridGrid::NakBsplineBoundaryCombigridGrid(size_t dim, size_t degree,
                                                                 level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SNakBsplineBoundaryCombigridBase(degree)),
      boundaryLevel(boundaryLevel) {}

NakBsplineBoundaryCombigridGrid::~NakBsplineBoundaryCombigridGrid() {}

sgpp::base::GridType NakBsplineBoundaryCombigridGrid::getType() {
  return sgpp::base::GridType::NakBsplineBoundaryCombigrid;
}

SBasis& NakBsplineBoundaryCombigridGrid::getBasis() { return *basis_; }

size_t NakBsplineBoundaryCombigridGrid::getDegree() { return this->degree; }

Grid* NakBsplineBoundaryCombigridGrid::unserialize(std::istream& istr) {
  return new NakBsplineBoundaryCombigridGrid(istr);
}

void NakBsplineBoundaryCombigridGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NakBsplineBoundaryCombigridGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
