// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINETRUNCATEDBOUNDARYGRID_HPP
#define BSPLINETRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Grid with Bspline basis functions with boundaries, pentagon cut
 */
class BsplineBoundaryGrid : public Grid {
 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   *
   * @param istr inputstream that contains the grid information
   */
  explicit BsplineBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with Bspline basis functions with boundaries, pentagon cut
   *
   * @param dim the dimension of the grid
   * @param degree the bspline's degree
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  BsplineBoundaryGrid(size_t dim,
                      size_t degree,
                      level_t boundaryLevel = 1);

  /**
   * Destructor.
   */
  ~BsplineBoundaryGrid() override;

  /**
   * @return string that identifies the grid type uniquely
   */
  SGPP::base::GridType getType() override;

  /**
   * @return B-spline basis
   */
  const SBasis& getBasis() override;

  /**
   * @return pointer to a GridGenerator object
   */
  GridGenerator& getGenerator() override;

  /**
   * reads a grid out of a string
   *
   * @param istr string that contains the grid information
   * @return grid
   */
  static std::unique_ptr<Grid> unserialize(std::istream& istr);

  /**
   * Serializes the grid.
   *
   * @param ostr stream to which the grid is written
   */
  void serialize(std::ostream& ostr) override;

  /**
   * @return B-spline degree
   */
  virtual size_t getDegree();

 protected:
  /// grid generator
  BoundaryGridGenerator generator;
  /// B-spline degree
  size_t degree;
  /// B-spline basis
  const SBsplineBoundaryBase* basis_;
  /// 1 + how much levels the boundary is coarser than the main axes
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* BSPLINETRUNCATEDBOUNDARYGRID_HPP */
