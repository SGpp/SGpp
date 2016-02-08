// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>


#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

LinearStretchedGrid::LinearStretchedGrid(std::istream& istr) :
  Grid(istr) {

}

LinearStretchedGrid::LinearStretchedGrid(size_t dim) :
  Grid(dim) {
}

LinearStretchedGrid::LinearStretchedGrid(Stretching& BB) :
  Grid(BB) {
}

LinearStretchedGrid::~LinearStretchedGrid() {
}

SGPP::base::GridType LinearStretchedGrid::getType() {
  return SGPP::base::GridType::LinearStretched;
}

const SBasis& LinearStretchedGrid::getBasis() {
  static SLinearStretchedBase basis;
  return basis;
}

Grid* LinearStretchedGrid::unserialize(std::istream& istr) {
  return new LinearStretchedGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearStretchedGrid::createGridGenerator() {
  return new StandardGridGenerator(this->storage);
}
}
}
