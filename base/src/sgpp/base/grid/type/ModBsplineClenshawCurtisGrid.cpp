// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

ModBsplineClenshawCurtisGrid::ModBsplineClenshawCurtisGrid(std::istream& istr) :
  Grid(istr),
  degree(1 << 16),
  basis_(NULL) {
  istr >> degree;
}


ModBsplineClenshawCurtisGrid::ModBsplineClenshawCurtisGrid(size_t dim,
    size_t degree) :
  Grid(dim),
  degree(degree),
  basis_(NULL) {
}

ModBsplineClenshawCurtisGrid::~ModBsplineClenshawCurtisGrid() {
  if (basis_ != NULL) {
    delete basis_;
  }
}

SGPP::base::GridType ModBsplineClenshawCurtisGrid::getType() {
  return SGPP::base::GridType::ModBsplineClenshawCurtis;
}

const SBasis& ModBsplineClenshawCurtisGrid::getBasis() {
  if (basis_ == NULL) {
    basis_ = new SBsplineModifiedClenshawCurtisBase(degree);
  }

  return *basis_;
}

size_t ModBsplineClenshawCurtisGrid::getDegree() {
  return this->degree;
}

Grid* ModBsplineClenshawCurtisGrid::unserialize(std::istream& istr) {
  return new ModBsplineClenshawCurtisGrid(istr);
}

void ModBsplineClenshawCurtisGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModBsplineClenshawCurtisGrid::createGridGenerator() {
  return new StandardGridGenerator(this->storage);
}

}  // namespace base
}  // namespace SGPP
