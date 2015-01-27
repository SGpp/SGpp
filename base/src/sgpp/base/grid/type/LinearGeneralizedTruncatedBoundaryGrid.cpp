// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>
#include "LinearGeneralizedTruncatedBoundaryGrid.hpp"
#include "../generation/GeneralizedTruncatedBoundaryGridGenerator.hpp"


namespace SGPP {
  namespace base {

    LinearGeneralizedTruncatedBoundaryGrid::LinearGeneralizedTruncatedBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearGeneralizedTruncatedBoundaryGrid::LinearGeneralizedTruncatedBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearGeneralizedTruncatedBoundaryGrid::LinearGeneralizedTruncatedBoundaryGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearGeneralizedTruncatedBoundaryGrid::~LinearGeneralizedTruncatedBoundaryGrid() {
    }

    const char* LinearGeneralizedTruncatedBoundaryGrid::getType() {
      return "linearGeneralizedTruncatedBoundary";
    }

    const SBasis& LinearGeneralizedTruncatedBoundaryGrid::getBasis() {
      static SLinearBoundaryBase basis;
      return basis;
    }

    Grid* LinearGeneralizedTruncatedBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearGeneralizedTruncatedBoundaryGrid(istr);
    }
    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearGeneralizedTruncatedBoundaryGrid::createGridGenerator() {
      return new GeneralizedTruncatedBoundaryGridGenerator(this->storage);
    }

  }
}
