// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODBSPLINECLENSHAWCURTISGRID_HPP
#define MODBSPLINECLENSHAWCURTISGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Grid with modified Clenshaw-Curtis Bspline basis functions
 */
class ModBsplineClenshawCurtisGrid : public Grid {
 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   *
   * @param istr inputstream that contains the grid information
   */
  explicit ModBsplineClenshawCurtisGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with modified Clenshaw-Curtis Bspline basis functions
   *
   * @param dim the dimension of the grid
     * @param degree the bspline's degree
   */
  explicit ModBsplineClenshawCurtisGrid(size_t dim, size_t degree);

  /**
   * Destructor.
   */
  ~ModBsplineClenshawCurtisGrid() override;

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
  std::unique_ptr<GridGenerator> createGridGenerator() override;

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
  /// B-spline degree
  size_t degree;
  /// B-spline basis
  const SBsplineModifiedClenshawCurtisBase* basis_;
};

}  // namespace base
}  // namespace SGPP

#endif /* MODBSPLINECLENSHAWCURTISGRID_HPP */
