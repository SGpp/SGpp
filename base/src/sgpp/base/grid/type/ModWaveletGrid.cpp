// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

ModWaveletGrid::ModWaveletGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

ModWaveletGrid::ModWaveletGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

ModWaveletGrid::~ModWaveletGrid() {
}

sgpp::base::GridType ModWaveletGrid::getType() {
  return sgpp::base::GridType::ModWavelet;
}

const SBasis& ModWaveletGrid::getBasis() {
  static SWaveletModifiedBase basis;
  return basis;
}

Grid* ModWaveletGrid::unserialize(std::istream& istr) {
  return new ModWaveletGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModWaveletGrid::getGenerator() {
  return generator;
}



}  // namespace base
}  // namespace sgpp
