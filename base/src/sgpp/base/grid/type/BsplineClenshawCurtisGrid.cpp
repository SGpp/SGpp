// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>

#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(std::istream& istr) :
      Grid(istr),
      degree(1 << 16),
      basis_(NULL),
      boundaryLevel(0) {
      istr >> degree;
      istr >> boundaryLevel;
    }

    BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(size_t dim,
        size_t degree,
        level_t boundaryLevel) :
      Grid(dim),
      degree(degree),
      basis_(NULL),
      boundaryLevel(boundaryLevel) {
    }

    BsplineClenshawCurtisGrid::~BsplineClenshawCurtisGrid() {
      if (basis_ != NULL) {
        delete basis_;
      }
    }

    SGPP::base::GridType BsplineClenshawCurtisGrid::getType() {
      return SGPP::base::GridType::BsplineClenshawCurtis;
    }

    const SBasis& BsplineClenshawCurtisGrid::getBasis() {
      if (basis_ == NULL) {
        basis_ = new SBsplineClenshawCurtisBase(degree);
      }

      return *basis_;
    }

    size_t BsplineClenshawCurtisGrid::getDegree() {
      return this->degree;
    }

    Grid* BsplineClenshawCurtisGrid::unserialize(std::istream& istr) {
      return new BsplineClenshawCurtisGrid(istr);
    }

    void BsplineClenshawCurtisGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
      ostr << boundaryLevel << std::endl;
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* BsplineClenshawCurtisGrid::createGridGenerator() {
      return new BoundaryGridGenerator(this->storage, boundaryLevel);
    }

  }
}
