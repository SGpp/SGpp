// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>


#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>


#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    LinearGrid::LinearGrid(std::istream& istr) :
      Grid(istr) {

    }

    LinearGrid::LinearGrid(size_t dim) :
      Grid(dim) {
    }

    LinearGrid::LinearGrid(BoundingBox& BB) :
      Grid(BB) {
    }

    LinearGrid::~LinearGrid() {
    }

    SGPP::base::GridType LinearGrid::getType() {
      return SGPP::base::GridType::Linear;
    }

    const SBasis& LinearGrid::getBasis() {
      static SLinearBase basis;
      return basis;
    }

    Grid* LinearGrid::unserialize(std::istream& istr) {
      return new LinearGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}
