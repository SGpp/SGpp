// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/LinearGridStencil.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>
#include <exception>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    LinearGridStencil::LinearGridStencil(std::istream& istr) : GridStencil(istr) {

    }

    LinearGridStencil::LinearGridStencil(size_t dim) : GridStencil(dim) {
      this->storage = new GridStorage(dim);
    }

    LinearGridStencil::LinearGridStencil(BoundingBox& BB) : GridStencil(BB) {
      this->storage = new GridStorage(BB);
    }

    LinearGridStencil::~LinearGridStencil() {
    }

    const char* LinearGridStencil::getType() {
      return "linearstencil";
    }

    const SBasis& LinearGridStencil::getBasis(){
		throw new factory_exception("Not implemented");
		// it should never get so far, code just for compilation reasons
		// If there will be a meaningful basis, this following lines should be changed
		static SLinearBase basis;
		return basis;
	}

    Grid* LinearGridStencil::unserialize(std::istream& istr) {
      return new LinearGridStencil(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearGridStencil::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}