// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>



#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    ModLinearGrid::ModLinearGrid(std::istream& istr) :
      Grid(istr) {
    }

    ModLinearGrid::ModLinearGrid(size_t dim) :
      Grid(dim) {
    }

    ModLinearGrid::~ModLinearGrid() {
    }

    SGPP::base::GridType ModLinearGrid::getType() {
      return SGPP::base::GridType::ModLinear;
    }

    const SBasis& ModLinearGrid::getBasis() {
      static SLinearModifiedBase basis;
      return basis;
    }

    Grid* ModLinearGrid::unserialize(std::istream& istr) {
      return new ModLinearGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* ModLinearGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}
