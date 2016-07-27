// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

WaveletGrid::WaveletGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

WaveletGrid::WaveletGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

WaveletGrid::~WaveletGrid() {
}

sgpp::base::GridType WaveletGrid::getType() {
  return sgpp::base::GridType::Wavelet;
}

SBasis& WaveletGrid::getBasis() {
  static SWaveletBase basis;
  return basis;
}

std::unique_ptr<Grid> WaveletGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new WaveletGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& WaveletGrid::getGenerator() {
  return generator;
}



}  // namespace base
}  // namespace sgpp
