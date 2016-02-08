// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODFUNDAMENTALSPLINEGRID_HPP
#define MODFUNDAMENTALSPLINEGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Grid with modified fundamental spline basis functions
 */
class ModFundamentalSplineGrid : public Grid {
 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   *
   * @param istr inputstream that contains the grid information
   */
  explicit ModFundamentalSplineGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with modified fundamental spline basis functions
   *
   * @param dim the dimension of the grid
   * @param degree fundamental spline degree
   */
  explicit ModFundamentalSplineGrid(size_t dim, size_t degree);

  /**
   * Destructor.
   */
  ~ModFundamentalSplineGrid() override;

  /**
   * @return string that identifies the grid type uniquely
   */
  SGPP::base::GridType getType() override;

  /**
   * @return fundamental spline basis
   */
  const SBasis& getBasis() override;

  /**
   * @return pointer to a GridGenerator object
   */
  GridGenerator* createGridGenerator() override;

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
   */
  void serialize(std::ostream& ostr) override;

  /**
   * @return B-spline degree
   */
  virtual size_t getDegree();

 protected:
  /// fundamental spline degree
  size_t degree;
  /// fundamental spline basis
  const SFundamentalSplineModifiedBase* basis_;
};

}  // namespace base
}  // namespace SGPP

#endif /* MODFUNDAMENTALSPLINEGRID_HPP */
