/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Gerrit Buse (buse@in.tum.de)

#include <sgpp/base/grid/type/GridStencil.hpp>

namespace sg {
  namespace base {

    GridStencil::GridStencil(std::istream& istr)
      : Grid(istr), surplusStencil(64), neighborStencil(64), weightStencil(64) {
    }


    GridStencil::GridStencil(size_t dim)
      : surplusStencil(64), neighborStencil(64), weightStencil(64) {
    }


    GridStencil::GridStencil(BoundingBox& BB)
      : surplusStencil(64), neighborStencil(64), weightStencil(64) {
    }


    GridStencil::~GridStencil() {
    }



    const GridStencil::IndexStencil&
    GridStencil::getSurplusStencil() const {
      return surplusStencil;
    }


    const GridStencil::IndexStencil&
    GridStencil::getNeighborStencil() const {
      return neighborStencil;
    }


    const GridStencil::WeightStencil&
    GridStencil::getWeightStencil() const {
      return weightStencil;
    }


  }
}

