// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>


namespace SGPP {
  namespace base {

    LinearBoundaryGrid::LinearBoundaryGrid(std::istream& istr) :
      Grid(istr),
      boundaryLevel(0) {
      istr >> boundaryLevel;
    }

    LinearBoundaryGrid::LinearBoundaryGrid(size_t dim, level_t boundaryLevel) :
      Grid(dim),
      boundaryLevel(boundaryLevel) {
    }

    LinearBoundaryGrid::LinearBoundaryGrid(BoundingBox& BB,
                                           level_t boundaryLevel) :
      Grid(BB),
      boundaryLevel(boundaryLevel) {
    }

    LinearBoundaryGrid::~LinearBoundaryGrid() {
    }

    SGPP::base::GridType LinearBoundaryGrid::getType() {
      return SGPP::base::GridType::LinearBoundary;
    }

    const SBasis& LinearBoundaryGrid::getBasis() {
      static SLinearBoundaryBase basis;
      return basis;
    }

    Grid* LinearBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearBoundaryGrid(istr);
    }

    void LinearBoundaryGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << boundaryLevel << std::endl;
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearBoundaryGrid::createGridGenerator() {
      return new BoundaryGridGenerator(this->storage, boundaryLevel);
    }


  }
}
