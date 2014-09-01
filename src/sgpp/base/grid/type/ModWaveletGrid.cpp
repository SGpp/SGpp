/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/ModWaveletGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg {
  namespace base {

    ModWaveletGrid::ModWaveletGrid(std::istream& istr) : Grid(istr) {
    }

    ModWaveletGrid::ModWaveletGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    ModWaveletGrid::~ModWaveletGrid() {
    }

    const char* ModWaveletGrid::getType() {
      return "modWavelet";
    }

    Grid* ModWaveletGrid::unserialize(std::istream& istr) {
      return new ModWaveletGrid(istr);
    }

    GridGenerator* ModWaveletGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }

  }
}
