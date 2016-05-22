// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STRETCHING_HPP
#define STRETCHING_HPP

#define LOOKUPSIZE 2047
#define LOOKUPMAX 11

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/globaldef.hpp>

#include <time.h>

#include <cstddef>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace base {

struct Stretching1D {
  std::string type;
  double x_0, xsi;
  // Index 0 : current position
  // Index 1: left neighbor position
  // Index 2: right neighbor position
  double lookup[LOOKUPSIZE][3];
};

/**
 *
 * Stretching can be done in four different ways:
 *
 * 1. No Stretching operation is done
 * 2. Stretching is done via log/exp calculation
 * 3. Sinh Stretching
 * 4. Discrete Stretching (No function used, grid points are directly taken)
 */

class Stretching : public BoundingBox {
 private:
  std::vector<Stretching1D> stretching1Ds;
  std::vector<int> discreteVectorLevel;
  std::string stretchingMode;

  /*
   * generates the lookup table up to level 11 for all dimensions, the rest is calculated via the
   * Stretching type given.
   */
  void generateLookupTable();

  /*
   * generates the positions of the left and right neighbors of one specific level and index point.
   * Written to the lookup table
   */
  void generateLeftRightArrays(Stretching1D& str1D, size_t d);

  /*
   * The main Stretching function uses the formula
   *
   * \f$ x_i = f^{-1}(f(a) + \frac{i}{n}.(f(b)-f(a)) ), i=0,...,n \f$
   * where f is the transformation function and [a,b] are the intervals the transformation is defined.
   *
   * Calculates the streched coordinate with respect to the given parameters
   *Note to self: may take other than double, take care
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   */
  double stretchingXform(int level, int index, size_t d) const;

  /*
   * calculates the logarithm transform in the given dimension
   *
   * \f$ f = log(x) \text{and} f^{-1}=e^x \f$
   *
   * @param d describes in which dimension the Stretching occurs.
   */
  void logXform(Stretching1D& str1D, size_t d) const;

  /*
   * overloaded
   * calculates the logarithm transform for an arbitrary input
   *
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   */
  double logXform(int level, int index, size_t d) const;

  /*
   * calculates the leentvar transform in the given dimension
   *
   * \f$ f= sinh_{-1}(\xi(x-x_0) \text{and} f^{-1}=frac{1}{\xi}sinh(y)+x_0\f$
   *
   * @param d describes in which dimension the Stretching occurs.
   */
  void leentvaarXform(Stretching1D& str1D, size_t d) const;

  /*
   * overloaded
   * calculates the leentvar transform in the given dimension
   *
   * \f$ f= sinh_{-1}(\xi(x-x_0) \text{and} f^{-1}=frac{1}{\xi}sinh(y)+x_0\f$
   *
   * @param d describes in which dimension the Stretching occurs.
   */
  double leentvaarXform(int level, int index, size_t d) const;

  /*
   * makes no Stretching on the dimension, the points are equidistant
   *
   * @param d describes in which dimension the Stretching occurs.
   */
  void noXform(Stretching1D& str1D, size_t d) const;

  /*
   * overloaded
   * returns the position of the node determined by the input arguments, on the non-streched grid.
   *
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   *
   */
  double noXform(int level, int index, size_t d) const;

  /*
   * calculates the lookup table index
   *
   * @param level level of the point
   * @param index index of the point
   */
  int calculateLookupIndex(int level, int index) const;

  /*
   * gets the discrete points that are stretched and creates a lookup table
   * @param vec double vector containing the discrete points, includes left and right boundary
   * @param stretch1d Stretch1D to write on the type of stretching and lookup table
   * @param d dimension of the node
   * @param discreteVectorLevel how many levels the grid vector points form
   */

  void parseVectorToLookupTable(std::vector<double>& vec,
                                Stretching1D& stretch1d,
                                size_t d, int& discreteVectorLevel) const;

  /*
   * calculates the left and right neighbor index and levels
   * @param level current level of the node
   * @param index current index of the node
   * @param leftLevel left neighbor level
   * @param leftIndex left neighbor index
   * @param rightLevel right neighbor level
   * @param rightIndex right neighbor index
   */
  void calculateNeighborSpecs(int level, int index, int& leftLevel,
                              int& leftIndex, int& rightLevel, int& rightIndex) const;

 public:
  /**
   * Constructor for Stretching
   *
   * initializes the Stretching using the boundaries with the input type array given
   * for each dimension
   *
   * @param boundaries BoundingBox1D struct to get the boundary values for Stretching
   * @param stretching1Ds array to define the stretching type
   */
  Stretching(const std::vector<BoundingBox1D>& boundaries,
             const std::vector<Stretching1D>& stretching1Ds);

  /**
   * Constructor for Stretching
   *
   * initializes the Stretching using the coordinates given. (For Janos' request)
   *
   * @param dimension number of the dimensions used with the grid
   * @param coordinates vector<double> array to get the boundaries, as well as the coordinates
   * of the specific level the vector defines for each dimension.
   */
  Stretching(size_t dimension, std::vector<double>* coordinates);

  /*
   * Gets the node position defined in the input parameters
   *
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   */
  double getCoordinate(int level, int index, size_t d) const;

  /*
   * Returns the type of the Stretching of the dimension given
   *
   * @param d dimension index of the type array
   *
   */
  Stretching1D getStretching1D(size_t d) const;

  /*
   * Prints the lookup table generated. Used for debugging
   */
  void printLookupTable() const;

  /*
   * Given the level, index and dimension, returns the lookup entries for left, current and right positions.
   *
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   * @param posc the point itself
   * @param posl left point
   * @param posr right point
   */
  void getAdjacentPositions(int level, int index, size_t d,
                            double& posc, double& posl, double& posr) const;

  /*
   * returns the stretchingMode, can be "analytic" or "discrete"
   */
  std::string getStretchingMode() const;

  /*
   * returns the discreteVector for the discrete grid structure.
   * used for serializing the grid
   */
  std::vector<double>* getDiscreteVector(bool bSort) const;

  /*
   * returns the discretevectorlevel array of size nDim
   */
  std::vector<int> getDiscreteVectorLevel() const;

  void calculateNeighborLookup(int maxlevel) const;
};

}  // namespace base
}  // namespace sgpp

#endif /* STRETCHING_HPP */
