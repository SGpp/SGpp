// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef NAKBSPLINETRUNCATEDBOUNDARYCOMBIGRIDGRID_HPP
#define NAKBSPLINETRUNCATEDBOUNDARYCOMBIGRIDGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryCombigridBasis.hpp>

namespace sgpp {
namespace base {

/**
 * Grid with B-spline basis functions with not-a-knot boundary conditions
 */
class NakBsplineBoundaryCombigridGrid : public Grid {
 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   *
   * @param istr inputstream that contains the grid information
   */
  explicit NakBsplineBoundaryCombigridGrid(std::istream& istr);

 public:
  /**
   * Constructor
   *
   * @param dim the dimension of the grid
   * @param degree the bspline's degree
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  NakBsplineBoundaryCombigridGrid(size_t dim, size_t degree, level_t boundaryLevel = 0);

  /**
   * Destructor.
   */
  ~NakBsplineBoundaryCombigridGrid() override;

  /**
   * @return string that identifies the grid type uniquely
   */
  sgpp::base::GridType getType() override;

  /**
   * @return B-spline basis
   */
  SBasis& getBasis() override;

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
  static Grid* unserialize(std::istream& istr);

  /**
   * Serializes the grid.
   *
   * @param ostr stream to which the grid is written
   * @param version the serialization version of the file
   *
   */
  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

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
  std::unique_ptr<SNakBsplineBoundaryCombigridBase> basis_;
  /// 1 + how much levels the boundary is coarser than the main axes
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace sgpp

#endif /* NAKBSPLINETRUNCATEDBOUNDARYCOMBIGRIDGRID_HPP */
