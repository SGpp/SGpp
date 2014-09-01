/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/LinearClenshawCurtisGrid.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include <iostream>

namespace sg {
  namespace base {

    LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream& istr) : Grid(istr) {
    }

    LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim,
        const CosineTable* cosine_table) :
      cosine_table(cosine_table) {
      this->storage = new GridStorage(dim);
    }

    LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid() {
    }

    const char* LinearClenshawCurtisGrid::getType() {
      return "linearClenshawCurtis";
    }

    base::Grid* LinearClenshawCurtisGrid::unserialize(std::istream& istr) {
      return new LinearClenshawCurtisGrid(istr);
    }

    base::GridGenerator* LinearClenshawCurtisGrid::createGridGenerator() {
      return new TrapezoidBoundaryGridGenerator(this->storage);
    }

    const CosineTable* LinearClenshawCurtisGrid::getCosineTable() const {
      return this->cosine_table;
    }

    void LinearClenshawCurtisGrid::setCosineTable(const CosineTable* cosine_table) {
      this->cosine_table = cosine_table;
    }

  }
}
