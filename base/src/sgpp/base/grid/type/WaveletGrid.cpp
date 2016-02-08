// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>



#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

WaveletGrid::WaveletGrid(std::istream& istr) :
  Grid(istr) {
}

WaveletGrid::WaveletGrid(size_t dim) :
  Grid(dim) {
}

WaveletGrid::~WaveletGrid() {
}

SGPP::base::GridType WaveletGrid::getType() {
  return SGPP::base::GridType::Wavelet;
}

const SBasis& WaveletGrid::getBasis() {
  static SWaveletBase basis;
  return basis;
}

Grid* WaveletGrid::unserialize(std::istream& istr) {
  return new WaveletGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* WaveletGrid::createGridGenerator() {
  return new StandardGridGenerator(this->storage);
}



}  // namespace base
}  // namespace SGPP
