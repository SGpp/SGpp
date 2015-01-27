// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>
#include "LinearStretchedTruncatedBoundaryGrid.hpp"
#include "../generation/StretchedTruncatedBoundaryGridGenerator.hpp"


namespace SGPP {
  namespace base {

    LinearStretchedTruncatedBoundaryGrid::LinearStretchedTruncatedBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearStretchedTruncatedBoundaryGrid::LinearStretchedTruncatedBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearStretchedTruncatedBoundaryGrid::LinearStretchedTruncatedBoundaryGrid(Stretching& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearStretchedTruncatedBoundaryGrid::~LinearStretchedTruncatedBoundaryGrid() {
    }

    const char* LinearStretchedTruncatedBoundaryGrid::getType() {
      return "linearStretchedTruncatedBoundary";
    }

    const SBasis& LinearStretchedTruncatedBoundaryGrid::getBasis() {
      static SLinearStretchedBoundaryBase basis;
      return basis;
    }

    Grid* LinearStretchedTruncatedBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearStretchedTruncatedBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearStretchedTruncatedBoundaryGrid::createGridGenerator() {
      return new StretchedTruncatedBoundaryGridGenerator(this->storage);
    }

  }
}