// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/WaveletTruncatedBoundaryGrid.hpp>

#include <sgpp/base/grid/generation/TruncatedBoundaryGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>



#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    WaveletTruncatedBoundaryGrid::WaveletTruncatedBoundaryGrid(std::istream& istr) : Grid(istr) {
    }

    WaveletTruncatedBoundaryGrid::WaveletTruncatedBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    WaveletTruncatedBoundaryGrid::~WaveletTruncatedBoundaryGrid() {
    }

    const char* WaveletTruncatedBoundaryGrid::getType() {
      return "waveletTruncatedBoundary";
    }

    const SBasis& WaveletTruncatedBoundaryGrid::getBasis() {
      static SWaveletBoundaryBase basis;
      return basis;
    }

    Grid* WaveletTruncatedBoundaryGrid::unserialize(std::istream& istr) {
      return new WaveletTruncatedBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* WaveletTruncatedBoundaryGrid::createGridGenerator() {
      return new TruncatedBoundaryGridGenerator(this->storage);
    }



  }
}
