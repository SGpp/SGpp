// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BOUNDINGBOX_HPP
#define BOUNDINGBOX_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/storage/hashmap/SerializationVersion.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * struct that defines the boundaries for one specific dimension
 */
struct BoundingBox1D {
  /// left boundary
  double leftBoundary;
  /// right boundary
  double rightBoundary;
  /// whether to use Dirichlet boundaries on the left boundary
  bool bDirichletLeft;
  /// whether to use Dirichlet boundaries on the right boundary
  bool bDirichletRight;

  /**
   * Default constructor initializing leftBoundary = 0, rightBoundary = 1, and
   * bDirichletLeft = bDirichletLeft = false.
   */
  BoundingBox1D() : BoundingBox1D(0.0, 1.0, false, false) {}

  /**
   * Constructor initializing bDirichletLeft = bDirichletLeft = false.
   *
   * @param leftBoundary  left boundary position
   * @param rightBoundary right boundary position
   */
  BoundingBox1D(double leftBoundary, double rightBoundary) :
    BoundingBox1D(leftBoundary, rightBoundary, false, false) {}

  /**
   * Constructor.
   *
   * @param leftBoundary    left boundary position
   * @param rightBoundary   right boundary position
   * @param bDirichletLeft  whether to use Dirichlet boundaries on the left boundary
   * @param bDirichletRight whether to use Dirichlet boundaries on the right boundary
   */
  BoundingBox1D(double leftBoundary, double rightBoundary,
                bool bDirichletLeft, bool bDirichletRight) :
                  leftBoundary(leftBoundary), rightBoundary(rightBoundary),
                  bDirichletLeft(bDirichletLeft), bDirichletRight(bDirichletRight) {}
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
  std::vector<BoundingBox1D> boundingBox1Ds;

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
   * @param boundingBox1Ds array that contains all boundaries
   */
  explicit BoundingBox(const std::vector<BoundingBox1D>& boundingBox1Ds);

  /**
   * Destructor.
   */
  virtual ~BoundingBox();

  /**
   * Sets left and right boundary for a specific dimension.
   *
   * @param d             the dimension in which the boundary should be changed
   * @param boundingBox1D reference to a BoundingBox1D object that contains the new boundaries
   */
  void setBoundary(size_t d, const BoundingBox1D& boundingBox1D);

  /**
   * Returns the left and right boundary for a specific dimension.
   *
   * @param d   the dimension in which the boundary should be read
   * @return    a BoundingBox1D object that contains the boundaries
   */
  inline const BoundingBox1D& getBoundary(size_t d) const {
    return boundingBox1Ds[d];
  }

  /**
   * Returns the number of dimensions of this bounding box.
   *
   * @return number of dimensions
   */
  size_t getDimension() const;

  /**
   * Calculates the width of the interval in one dimension.
   *
   * @param d the dimension in which the width of the interval should be determined
   * @return  width of the interval
   */
  inline double getIntervalWidth(size_t d) const {
    return boundingBox1Ds[d].rightBoundary - boundingBox1Ds[d].leftBoundary;
  }

  /**
   * Returns the offset in positive x-direction of the interval in one dimension.
   *
   * @param d dimension in which the offset of the interval should be determined
   * @return  offset in positive x-direction of the interval
   */
  inline double getIntervalOffset(size_t d) const {
    return boundingBox1Ds[d].leftBoundary;
  }

  /**
   * Determine if this bounding box describes the unit cube \f$[0, 1]^d\f$.
   *
   * @return true if this bounding box is the unit cube, otherwise false
   */
  bool isUnitCube() const;

  /**
   * Transform a point in the unit cube \f$[0, 1]^d\f$ to a point in the BoundingBox.
   *
   * @param[in,out] point point to be transformed in-place
   */
  inline void transformPointToBoundingBox(DataVector& point) const {
    for (size_t d = 0; d < dimension; d++) {
      point[d] = transformPointToBoundingBox(d, point[d]);
    }
  }

  /**
   * Transform a point in the unit interval \f$[0, 1]\f$ to a point in the BoundingBox in 1D.
   *
   * @param d       dimension
   * @param point   1D point in unit interval
   * @return        transformed 1D point in the BoundingBox
   */
  inline double transformPointToBoundingBox(size_t d, double point) const {
    return getIntervalOffset(d) + getIntervalWidth(d) * point;
  }

  /**
   * Transform a point in the BoundingBox to a point in the unit cube \f$[0, 1]^d\f$.
   *
   * @param[in,out] point point to be transformed in-place
   */
  inline void transformPointToUnitCube(DataVector& point) const {
    for (size_t d = 0; d < dimension; d++) {
      point[d] = transformPointToUnitCube(d, point[d]);
    }
  }

  /**
   * Transform a point in the BoundingBox to a point in the unit interval \f$[0, 1]\f$ in 1D.
   *
   * @param d       dimension
   * @param point   1D point in the BoundingBox
   * @return        transformed 1D point in the unit interval
   */
  inline double transformPointToUnitCube(size_t d, double point) const {
    return (point - getIntervalOffset(d)) / getIntervalWidth(d);
  }

  /**
   * Check whether the BoundingBox contains a given point.
   *
   * @param point   point to be checked
   * @return        whether the point is contained in the BoundingBox
   */
  inline bool isContainingPoint(DataVector& point) const {
    for (size_t d = 0; d < dimension; d++) {
      if (!isContainingPoint(d, point[d])) {
        return false;
      }
    }

    return true;
  }

  /**
   * Check whether the BoundingBox contains a given point in a specific dimension.
   *
   * @param d       dimension to be checked
   * @param point   1D point to be checked
   * @return        whether the point is contained in the BoundingBox in the given dimension
   */
  inline bool isContainingPoint(size_t d, double point) const {
    return (boundingBox1Ds[d].leftBoundary <= point) && (point <= boundingBox1Ds[d].rightBoundary);
  }

  /**
   * Determines if the interval in the specified dimension has left Dirichlet boundary conditions.
   *
   * @param d   the dimension for which the left boundary condition should be determined
   * @return    true if Dirichlet Boundary conditions, otherwise false
   */
  bool hasDirichletBoundaryLeft(size_t d) const;

  /**
   * Determines if the interval in the specified dimension has right Dirichlet boundary conditions.
   *
   * @param d   the dimension for which the right boundary condition should be determined
   * @return    true if Dirichlet Boundary conditions, otherwise false
   */
  bool hasDirichletBoundaryRight(size_t d) const;

  /**
   * Serialize the BoundingBox into a string.
   *
   * @param version the serialization version of the file
   * @return string that contains all BoundingBox information
   */
  virtual std::string serialize(int version = SERIALIZATION_VERSION) const;

  /**
   * Serialize the BoundingBox into a stream.
   *
   * @param ostream reference to a stream into that all BoundingBox information is written
   * @param version the serialization version of the file
   */
  virtual void serialize(std::ostream& ostream, int version = SERIALIZATION_VERSION) const;

  /**
   * Unserialize from a string.
   *
   * @param istr string which contains a serialized BoundingBox
   * @param version the serialization version of the file
   */
  virtual void unserialize(const std::string& istr, int version);

  /**
   * Unserialize from a stream.
   *
   * @param istr stream which contains a serialized BoundingBox
   * @param version the serialization version of the file
   */
  virtual void unserialize(std::istream& istr, int version);

  /**
   * Converts the BoundingBox to a string.
   *
   * @param text string to which the data is written
   */
  void toString(std::string& text) const;

  /**
   * Converts the BoundingBox to a string.
   *
   * @return string to which the data is written
   */
  std::string toString() const;
};

}  // namespace base
}  // namespace sgpp

#endif /* BOUNDINGBOX_HPP */
