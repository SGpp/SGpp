// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include "PolyTruncatedBoundaryGrid.hpp"

#include <sgpp/base/grid/generation/TruncatedBoundaryGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <iostream>

namespace SGPP {
  namespace base {

    PolyTruncatedBoundaryGrid::PolyTruncatedBoundaryGrid(std::istream& istr) :
      Grid(istr), degree(1 << 16), basis_(NULL) {
      istr >> degree;
    }

    PolyTruncatedBoundaryGrid::PolyTruncatedBoundaryGrid(size_t dim,
        size_t degree) :
      degree(degree), basis_(NULL) {
      this->storage = new GridStorage(dim);
    }

    PolyTruncatedBoundaryGrid::~PolyTruncatedBoundaryGrid() {
      if (basis_ != NULL) {
        delete basis_;
      }
    }

    const SBasis& PolyTruncatedBoundaryGrid::getBasis() {
      if (basis_ == NULL) {
        basis_ = new SPolyBoundaryBase(degree);
      }

      return *basis_;
    }

    const char* PolyTruncatedBoundaryGrid::getType() {
      return "polyTruncatedBoundary";
    }

    size_t PolyTruncatedBoundaryGrid::getDegree() const {
      return this->degree;
    }

    Grid* PolyTruncatedBoundaryGrid::unserialize(std::istream& istr) {
      return new PolyTruncatedBoundaryGrid(istr);
    }

    void PolyTruncatedBoundaryGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* PolyTruncatedBoundaryGrid::createGridGenerator() {
      return new TruncatedBoundaryGridGenerator(this->storage);
    }

  }
}
