// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELETGRID_HPP
#define WAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with wavelet base functions
 */
class WaveletGrid : public Grid {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit WaveletGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with wavelet base functions
   *
   * @param dim the dimension of the grid
   */
  explicit WaveletGrid(size_t dim);

  /**
   * Destructor
   */
  ~WaveletGrid() override;

  sgpp::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* WAVELETGRID_HPP */
