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
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  explicit WaveletBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~WaveletBoundaryGrid() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator* createGridGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);

 protected:
  /// 1 + how much levels the boundary is coarser than the main axes
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* WAVELETTRUNCATEDBOUNDARYGRID_HPP */
