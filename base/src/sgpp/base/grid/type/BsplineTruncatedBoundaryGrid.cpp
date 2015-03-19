// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/BsplineTruncatedBoundaryGrid.hpp>

#include <sgpp/base/grid/generation/TruncatedBoundaryGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    BsplineTruncatedBoundaryGrid::BsplineTruncatedBoundaryGrid(std::istream& istr) : Grid(istr), degree(1 << 16), basis_(NULL) {
      istr >> degree;
    }


    BsplineTruncatedBoundaryGrid::BsplineTruncatedBoundaryGrid(size_t dim, size_t degree) : degree(degree), basis_(NULL) {
      this->storage = new GridStorage(dim);
    }

    BsplineTruncatedBoundaryGrid::~BsplineTruncatedBoundaryGrid() {
      if (basis_ != NULL) {
        delete basis_;
      }
    }

    const char* BsplineTruncatedBoundaryGrid::getType() {
      return "bsplineTruncatedBoundary";
    }

    const SBasis& BsplineTruncatedBoundaryGrid::getBasis() {
      if (basis_ == NULL) {
        basis_ = new SBsplineBoundaryBase(degree);
      }

      return *basis_;
    }

    size_t BsplineTruncatedBoundaryGrid::getDegree() {
      return this->degree;
    }

    Grid* BsplineTruncatedBoundaryGrid::unserialize(std::istream& istr) {
      return new BsplineTruncatedBoundaryGrid(istr);
    }

    void BsplineTruncatedBoundaryGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* BsplineTruncatedBoundaryGrid::createGridGenerator() {
      return new TruncatedBoundaryGridGenerator(this->storage);
    }

  }
}
