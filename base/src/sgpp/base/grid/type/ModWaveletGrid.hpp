// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODWAVELETGRID_HPP
#define MODWAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with modified wavelet base functions
 */
class ModWaveletGrid : public Grid {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit ModWaveletGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with modified wavelet base functions
   *
   * @param dim the dimension of the grid
   */
  explicit ModWaveletGrid(size_t dim);

  /**
   * Destructor
   */
  ~ModWaveletGrid() override;

  sgpp::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* MODWAVELETGRID_HPP */
