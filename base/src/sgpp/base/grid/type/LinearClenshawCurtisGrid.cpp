// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>


namespace SGPP {
  namespace base {

    LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream& istr) : Grid(istr) {

    }

    LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid() {
    }

    const char* LinearClenshawCurtisGrid::getType() {
      return "linearClenshawCurtis";
    }

    const SBasis& LinearClenshawCurtisGrid::getBasis() {
      static SLinearClenshawCurtisBase basis;
      return basis;
    }

    Grid* LinearClenshawCurtisGrid::unserialize(std::istream& istr) {
      return new LinearClenshawCurtisGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearClenshawCurtisGrid::createGridGenerator() {
      return new BoundaryGridGenerator(this->storage);
    }


  }
}
