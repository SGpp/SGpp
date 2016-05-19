// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POLYTRUNCATEDBOUNDARYGRID_HPP
#define POLYTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * trapezoid boundary grid with polynomial base functions
 */
class PolyBoundaryGrid : public Grid {
 protected:
  explicit PolyBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with polynomial base functions
   *
   * @param dim the dimension of the grid
   * @param degree the max. polynom's degree
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  PolyBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~PolyBoundaryGrid() override;

  const SBasis& getBasis() override;
  sgpp::base::GridType getType() override;
  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
  size_t getDegree() const;

 protected:
  /// grid generator
  BoundaryGridGenerator generator;
  /// max. polynom's degree
  size_t degree;
  /// polynomial basis
  std::unique_ptr<SPolyBoundaryBase> basis_;
  /// 1 + how much levels the boundary is coarser than the main axes
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace sgpp

#endif /* POLYTRUNCATEDBOUNDARYGRID_HPP */
