/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Ricard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/PrewaveletGrid.hpp"

#include "base/grid/generation/PrewaveletGridGenerator.hpp"

#include "base/exception/factory_exception.hpp"

#include "base/basis/prewavelet/PrewaveletBasis.hpp"


namespace sg {
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
