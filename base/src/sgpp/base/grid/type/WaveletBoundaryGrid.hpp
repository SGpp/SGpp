// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELETTRUNCATEDBOUNDARYGRID_HPP
#define WAVELETTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with wavelet base functions with boundaries, pentagon cut
 */
class WaveletBoundaryGrid : public Grid {
 protected:
  explicit WaveletBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with wavelet base functions with boundaries, pentagon cut
   *
   * @param dim the dimension of the grid
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted (default = 1: boundary has same level
   *                      as main axes)
   */
  explicit WaveletBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~WaveletBoundaryGrid() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);

 protected:
  /// level at which the boundary points should be inserted
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* WAVELETTRUNCATEDBOUNDARYGRID_HPP */
