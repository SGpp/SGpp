/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"

namespace sg {
  namespace base {

    BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(std::istream& istr) :
      Grid(istr),
      degree(1 << 16) {
      istr >> degree;
    }

    BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(
      size_t dim, size_t degree, const CosineTable* cosine_table) :
      degree(degree),
      cosine_table(cosine_table) {
      this->storage = new GridStorage(dim);
    }

    BsplineClenshawCurtisGrid::~BsplineClenshawCurtisGrid() {
    }

    const char* BsplineClenshawCurtisGrid::getType() {
      return "BsplineClenshawCurtis";
    }

    size_t BsplineClenshawCurtisGrid::getDegree() {
      return this->degree;
    }

    const CosineTable* BsplineClenshawCurtisGrid::getCosineTable() const {
      return this->cosine_table;
    }

    void BsplineClenshawCurtisGrid::setCosineTable(const CosineTable* cosine_table) {
      this->cosine_table = cosine_table;
    }

    Grid* BsplineClenshawCurtisGrid::unserialize(std::istream& istr) {
      return new BsplineClenshawCurtisGrid(istr);
    }

    void BsplineClenshawCurtisGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }

    GridGenerator* BsplineClenshawCurtisGrid::createGridGenerator() {
      return new TrapezoidBoundaryGridGenerator(this->storage);
    }

  }
}
