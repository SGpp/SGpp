// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>

#include <sgpp/base/grid/generation/PrewaveletGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/basis/prewavelet/PrewaveletBasis.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    PrewaveletGrid::PrewaveletGrid(std::istream& istr) : Grid(istr) {
    }

    PrewaveletGrid::PrewaveletGrid(size_t dim) {
      this->storage = new GridStorage(dim);
      this->shadowStorage = new GridStorage(dim);
    }

    PrewaveletGrid::~PrewaveletGrid() {
    }

    const char* PrewaveletGrid::getType() {
      return "prewavelet";
    }

    const SBasis& PrewaveletGrid::getBasis(){
		static SPrewaveletBase basis;
		return basis;
	}

    Grid* PrewaveletGrid::unserialize(std::istream& istr) {
      return new PrewaveletGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* PrewaveletGrid::createGridGenerator() {
      return new PrewaveletGridGenerator(this->storage, this->shadowStorage);
    }


    GridStorage* PrewaveletGrid::getShadowStorage() {
      return this->shadowStorage;
    }

  }
}