/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/type/ModLinearGridStencil.hpp"

#include "base/grid/generation/StandardGridGenerator.hpp"


#include "base/exception/factory_exception.hpp"


#include <iostream>

namespace sg {
  namespace base {

    ModLinearGridStencil::ModLinearGridStencil(std::istream& istr) : GridStencil(istr) {

    }

    ModLinearGridStencil::ModLinearGridStencil(size_t dim) : GridStencil(dim) {
      this->storage = new GridStorage(dim);
    }

    ModLinearGridStencil::ModLinearGridStencil(BoundingBox& BB) : GridStencil(BB) {
      this->storage = new GridStorage(BB);
    }

    ModLinearGridStencil::~ModLinearGridStencil() {
    }

    const char* ModLinearGridStencil::getType() {
      return "modlinearstencil";
    }

    Grid* ModLinearGridStencil::unserialize(std::istream& istr) {
      return new ModLinearGridStencil(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* ModLinearGridStencil::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}
