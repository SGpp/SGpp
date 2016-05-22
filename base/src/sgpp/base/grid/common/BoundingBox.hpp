// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BOUNDINGBOX_HPP
#define BOUNDINGBOX_HPP

#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <string>


namespace sgpp {
namespace base {

/**
 * struct that defines the boundaries for one specific dimension
 */
struct BoundingBox1D {
  /// the left boundary
  double leftBoundary;
  /// the right boundary
  double rightBoundary;
  /// Use Dirichlet-Boundaries on the left boundary
  bool bDirichletLeft;
  /// Use Dirichlet-Boundaries on the right boundary
  bool bDirichletRight;
};

/**
 * This class implements the boundaries of the sparse grid.
 * Internally the grid is set up on a trivial cube.
 *
 * This class gives the class gives the opportunity to stretch
 * this cube in every dimension separately.
 */
class BoundingBox {
 protected:
  /// the number of dimensions used with the grid
  size_t dimension;
  /// Array that contains all left boundaries for all dimensions
  BoundingBox1D* dimensionBoundaries;
  /// is true if the complete Bounding Box desribes a trivial cube
  bool bTrivialCube;

 public:
  /**
   * Constructor for BoundingBox
   *
   * initializes the Bounding with a N-d trivial cube
   *
   * @param dimension number of the dimensions used with the grid
   */
  explicit BoundingBox(size_t dimension);

  /**
   * Constructor for BoundingBox
   *
   * initializes the Bounding with specific values for all dimensions
   *
   * @param dimension number of the dimensions used with the grid
   * @param boundaries array that contains all boundaries
   */
  BoundingBox(size_t dimension, const BoundingBox1D* boundaries);

  /**
   * Copy-Constructor
   *
   * initializes the Bounding with values of another bounding Box
   *
   * @param copyBoundingBox reference to a BoundingBox Object whose values are copied
   */
  BoundingBox(const BoundingBox& copyBoundingBox);

  /**
   * Desctructor
   */
  virtual ~BoundingBox();

  /**
   * set the left and right boundary for a specific dimension
   *
   * @param d the dimension in which the boundary should be changed
   * @param newBoundaries reference to a DimensionBoundary object, that contains the new boundaries
   */
  void setBoundary(size_t d, const BoundingBox1D& newBoundaries);

  /**
   * gets the left and right boundary for a specific dimension
   *
   * @param d the dimension in which the boundary should be read
   * @return a DimensionBoundary object, that contains the boundaries
   */
  BoundingBox1D getBoundary(size_t d) const;

  /**
   * gets the dimensions of the cube stored in this bounding box
   *
   * @return the number of dimensions
   */
  size_t getDimensions() const;

  /**
   * gets the width of the interval in one dimension
   *
   * @param d the dimension in which the width of the interval should be determined
   * @return the width of the interval
   */
  double getIntervalWidth(size_t d) const;

  /**
   * gets the offset in positive x-direction of the interval in one dimension
   *
   * @param d the dimension in which the offset of the interval should be determined
   * @return the offset in positive x-direction of the interval
   */
  double getIntervalOffset(size_t d) const;

  /**
   * Use this function to determine if this bounding box describes a trivial cube [0;1]^d
   *
   * @return true if this bounding boy is a trivial cube otherwise false
   */
  bool isTrivialCube() const;

  /**
   * Determines, if the interval in the specified dimension has left dirichlet boundary conditions
   *
   * @param dimension the dimension for which the left boundary condition should be determined
   * @return true if Dirichlet Boundary conditions, otherwise false
   */
  bool hasDirichletBoundaryLeft(size_t d) const;

  /**
   * Determines, if the interval in the specified dimension has right dirichlet boundary conditions
   *
   * @param dimension the dimension for which the right boundary condition should be determined
   * @return true if Dirichlet Boundary conditions, otherwise false
   */
  bool hasDirichletBoundaryRight(size_t d) const;

  void toString(std::string& text) const;

  std::string toString() const;
};

}  // namespace base
}  // namespace sgpp

#endif /* BOUNDINGBOX_HPP */
