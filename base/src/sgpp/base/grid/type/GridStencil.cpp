// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/GridStencil.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

GridStencil::GridStencil(std::istream& istr)
  : Grid(istr), surplusStencil(64), neighborStencil(64), weightStencil(64) {
}


GridStencil::GridStencil(size_t dim)
  : Grid(dim), surplusStencil(64), neighborStencil(64), weightStencil(64) {
}


GridStencil::GridStencil(BoundingBox& BB)
  : Grid(BB), surplusStencil(64), neighborStencil(64), weightStencil(64) {
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


}  // namespace base
}  // namespace sgpp
