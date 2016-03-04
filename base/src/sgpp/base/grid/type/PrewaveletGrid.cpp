// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/PrewaveletBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

PrewaveletGrid::PrewaveletGrid(std::istream& istr) :
  Grid(istr),
  shadowStorage(storage.getDimension()),
  generator(storage, shadowStorage) {
}

PrewaveletGrid::PrewaveletGrid(size_t dim) :
  Grid(dim),
  shadowStorage(dim),
  generator(storage, shadowStorage) {
}

PrewaveletGrid::~PrewaveletGrid() {
}

sgpp::base::GridType PrewaveletGrid::getType() {
  return base::GridType::Prewavelet;
}

const SBasis& PrewaveletGrid::getBasis() {
  static SPrewaveletBase basis;
  return basis;
}

std::unique_ptr<Grid> PrewaveletGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new PrewaveletGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& PrewaveletGrid::getGenerator() {
  return generator;
}


GridStorage& PrewaveletGrid::getShadowStorage() {
  return shadowStorage;
}

}  // namespace base
}  // namespace sgpp
