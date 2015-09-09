// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>
#include "LinearStretchedBoundaryGrid.hpp"
#include "../generation/StretchedBoundaryGridGenerator.hpp"


namespace SGPP {
  namespace base {

    LinearStretchedBoundaryGrid::LinearStretchedBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearStretchedBoundaryGrid::LinearStretchedBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearStretchedBoundaryGrid::LinearStretchedBoundaryGrid(Stretching& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearStretchedBoundaryGrid::~LinearStretchedBoundaryGrid() {
    }

    SGPP::base::GridType LinearStretchedBoundaryGrid::getType() {
      return SGPP::base::GridType::LinearStretchedBoundary;
    }

    const SBasis& LinearStretchedBoundaryGrid::getBasis() {
      static SLinearStretchedBoundaryBase basis;
      return basis;
    }

    Grid* LinearStretchedBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearStretchedBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearStretchedBoundaryGrid::createGridGenerator() {
      return new StretchedBoundaryGridGenerator(this->storage);
    }

  }
}
