// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


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

std::unique_ptr<Grid> ModBsplineClenshawCurtisGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new ModBsplineClenshawCurtisGrid(istr));
}

void ModBsplineClenshawCurtisGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
std::unique_ptr<GridGenerator> ModBsplineClenshawCurtisGrid::createGridGenerator() {
  return std::unique_ptr<GridGenerator>(new StandardGridGenerator(*this->storage));
}

}  // namespace base
}  // namespace SGPP
