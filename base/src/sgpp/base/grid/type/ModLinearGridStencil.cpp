/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>


#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/basis/modlinear/ModifiedLinearBasis.hpp>


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

    const SBasis& ModLinearGridStencil::getBasis(){
		throw new factory_exception("Not implemented");
		// it should never get so far, code just for compilation reasons
		// If there will be a meaningful basis, this following lines should be changed
		static SModLinearBase basis;
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
