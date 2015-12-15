// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    BsplineGrid::BsplineGrid(std::istream& istr) :
      Grid(istr),
      degree(1 << 16),
      basis_(NULL) {
      istr >> degree;
    }


    BsplineGrid::BsplineGrid(size_t dim, size_t degree) :
      Grid(dim),
      degree(degree),
      basis_(NULL) {
    }

    BsplineGrid::~BsplineGrid() {
      if (basis_ != NULL) {
        delete basis_;
      }
    }

    SGPP::base::GridType BsplineGrid::getType() {
      return SGPP::base::GridType::Bspline;
    }

    const SBasis& BsplineGrid::getBasis() {
      if (basis_ == NULL) {
        basis_ = new SBsplineBase(degree);
      }

      return *basis_;
    }

    size_t BsplineGrid::getDegree() {
      return this->degree;
    }

    Grid* BsplineGrid::unserialize(std::istream& istr) {
      return new BsplineGrid(istr);
    }

    void BsplineGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* BsplineGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }

  }
}
