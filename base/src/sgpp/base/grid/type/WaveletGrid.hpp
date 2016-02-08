// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELETGRID_HPP
#define WAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with wavelet base functions
 */
class WaveletGrid : public Grid {
 protected:
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

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace SGPP

#endif /* WAVELETGRID_HPP */
