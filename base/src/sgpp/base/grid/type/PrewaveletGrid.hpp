// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PREWAVELETGRID_HPP
#define PREWAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/PrewaveletGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with prewavelet base functions
 */
class PrewaveletGrid : public Grid {
 protected:
  GridStorage shadowStorage;
  /// grid generator
  PrewaveletGridGenerator generator;
  explicit PrewaveletGrid(std::istream& istr);

 public:
  explicit PrewaveletGrid(size_t dim);
  ~PrewaveletGrid() override;

  sgpp::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);

  /**
   * gets a reference to the GridStorage object
   *
   * @return reference to the GridStorage object
   */
  GridStorage& getShadowStorage();
};

}  // namespace base
}  // namespace sgpp

#endif /* PREWAVELETGRID_HPP */
