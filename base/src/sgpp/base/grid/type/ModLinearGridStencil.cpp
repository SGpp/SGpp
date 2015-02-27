// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>


#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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

    const SBasis& ModLinearGridStencil::getBasis() {
      throw new factory_exception("Not implemented");
      // it should never get so far, code just for compilation reasons
      // If there will be a meaningful basis, this following lines should be changed
      static SLinearModifiedBase basis;
      return basis;
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