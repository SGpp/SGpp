// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STRETCHING_HPP
#define STRETCHING_HPP

#define LOOKUPSIZE 2047
#define LOOKUPMAX 11

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>
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

  /**
   * Default constructor initializing type = "id", x_0 = 0.0, and xsi = 1.0.
   */
  Stretching1D() : Stretching1D("id") {}

  /**
   * Constructor initializing x_0 = 0.0, and xsi = 1.0.
   */
  explicit Stretching1D(const std::string& type) : type(type), x_0(0.0), xsi(1.0) {}
};

/**
 *
 * Stretching can be done in different ways:
 *
 * 1. No Stretching operation is done
 * 2. Clenshaw-Curtis Stretching
 * 3. Stretching is done via log/exp calculation
 * 4. Sinh Stretching
 * 5. Fitob Stretching
 * 6. Discrete Stretching (No function used, grid points are directly taken)
 */

class Stretching : public BoundingBox {
 private:
  std::vector<Stretching1D> stretching1Ds;
  std::vector<level_t> discreteVectorLevel;
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
   * Note to self: may take other than double, take care
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   */
  inline double stretchingXform(level_t level, index_t index, size_t d) const {
    const std::string& type = stretching1Ds[d].type;

    // refer to the lookup table
    if (level <= LOOKUPMAX) {
      return stretching1Ds[d].lookup[calculateLookupIndex(level, index)][0];
    } else if (type == "cc") {
      // calculate using the appropriate function
      return clenshawCurtisXform(level, index, d);
    } else if (type == "log") {
      return logXform(level, index, d);
    } else if (type == "sinh") {
      return leentvaarXform(level, index, d);
    } else if (type == "fitob") {
      // level - LOOKUPMAX is guaranteed to be > 0
      index_t pow2deltaL = static_cast<index_t>(1) << (level - LOOKUPMAX);
      double dIndex = static_cast<double>(index) / static_cast<double>(pow2deltaL);
      level_t leftLevel = LOOKUPMAX;
      level_t rightLevel = LOOKUPMAX;
      int leftIndex = static_cast<index_t>(std::floor(dIndex));
      int rightIndex = static_cast<index_t>(std::ceil(dIndex));
      bool leftContinue = (leftIndex % 2 == 0);
      bool rightContinue = (rightIndex % 2 == 0);

      while (leftContinue || rightContinue) {
        if (leftContinue) {
          leftIndex /= 2;
          leftLevel--;

          // if left index odd, or left level 0, stop calculating left part
          if (((leftIndex % 2 != 0)) || (leftLevel == 0)) {
            if (leftLevel == 0) {
              leftIndex = 0;
            }

            leftContinue = false;
          }
        }

        if (rightContinue) {
          rightIndex /= 2;
          rightLevel--;

          // if right index odd, or right level 0, stop calculating right part
          if (((rightIndex % 2) != 0) || (rightLevel == 0)) {
            if (rightLevel == 0) {
              rightIndex = 1;
            }

            rightContinue = false;
          }
        }
      }

      const double posl = getCoordinate(leftLevel, leftIndex, d);
      const double posr = getCoordinate(rightLevel, rightIndex, d);
      const double step = (posr - posl) / static_cast<double>(pow2deltaL);

      return posl + step * (dIndex - std::floor(dIndex)) * pow2deltaL;
    } else {
      return noXform(level, index, d);
    }
  }

  /*
   * calculates the Clenshaw-Curtis points in the given dimension
   *
   * @param d describes in which dimension the Stretching occurs
   */
  void clenshawCurtisXform(Stretching1D& str1D, size_t d) const;

  /*
   * overloaded
   * calculates the logarithm transform for an arbitrary input
   *
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   */
  inline double clenshawCurtisXform(level_t level, index_t index, size_t d) const {
    const double bbLeft = boundingBox1Ds[d].leftBoundary;
    const double bbWidth = boundingBox1Ds[d].rightBoundary - bbLeft;

    return bbLeft + bbWidth * ClenshawCurtisTable::getInstance().getPoint(level, index);
  }

  /*
   * calculates the logarithm transform in the given dimension
   *
   * \f$ f = log(x) \text{and} f^{-1}=e^x \f$
   *
   * @param d describes in which dimension the Stretching occurs
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
  inline double logXform(level_t level, index_t index, size_t d) const {
    const double bbLeftTransformed = std::log(boundingBox1Ds[d].leftBoundary);
    const double bbWidthTransformed = std::log(boundingBox1Ds[d].rightBoundary) - bbLeftTransformed;
    const index_t hInv = static_cast<index_t>(1) << level;

    return std::exp(bbLeftTransformed + bbWidthTransformed *
                    static_cast<double>(index) / static_cast<double>(hInv));
  }

  /*
   * calculates the leentvar transform in the given dimension
   *
   * \f$ f= sinh_{-1}(\xi(x-x_0) \text{and} f^{-1}=frac{1}{\xi}sinh(y)+x_0\f$
   *
   * @param d describes in which dimension the Stretching occurs
   */
  void leentvaarXform(Stretching1D& str1D, size_t d) const;

  /*
   * overloaded
   * calculates the leentvar transform in the given dimension
   *
   * \f$ f= sinh_{-1}(\xi(x-x_0) \text{and} f^{-1}=frac{1}{\xi}sinh(y)+x_0\f$
   *
   * @param d describes in which dimension the Stretching occurs
   */
  inline double leentvaarXform(level_t level, index_t index, size_t d) const {
    const double sinhArgumentLeft = ((boundingBox1Ds[d].leftBoundary -
        stretching1Ds[d].x_0) * stretching1Ds[d].xsi);
    const double sinhArgumentRight = ((boundingBox1Ds[d].rightBoundary -
        stretching1Ds[d].x_0) * stretching1Ds[d].xsi);
    const double bbLeftTransformed =
        std::log(sinhArgumentLeft + std::sqrt(sinhArgumentLeft * sinhArgumentLeft + 1));
    const double bbWidthTransformed =
        std::log(sinhArgumentRight + std::sqrt(sinhArgumentRight * sinhArgumentRight + 1)) -
        bbLeftTransformed;
    const index_t hInv = static_cast<index_t>(1) << level;

    return stretching1Ds[d].x_0 + 1.0 / stretching1Ds[d].xsi * std::sinh((
        bbLeftTransformed + bbWidthTransformed *
        static_cast<double>(index) / static_cast<double>(hInv)));
  }

  /*
   * makes no Stretching on the dimension, the points are equidistant
   *
   * @param d describes in which dimension the Stretching occurs
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
  inline double noXform(level_t level, index_t index, size_t d) const {
    const double bbLeft = boundingBox1Ds[d].leftBoundary;
    const double bbWidth = boundingBox1Ds[d].rightBoundary - bbLeft;
    const index_t hInv = static_cast<index_t>(1) << level;

    return bbLeft + bbWidth * static_cast<double>(index) / static_cast<double>(hInv);
  }

  /*
   * calculates the lookup table index
   *
   * @param level level of the point
   * @param index index of the point
   */
  size_t calculateLookupIndex(level_t level, index_t index) const;

  /*
   * gets the discrete points that are stretched and creates a lookup table
   * @param vec double vector containing the discrete points, includes left and right boundary
   * @param stretch1d Stretch1D to write on the type of stretching and lookup table
   * @param d dimension of the node
   * @param discreteVectorLevel how many levels the grid vector points form
   */

  void parseVectorToLookupTable(std::vector<double>& vec,
                                Stretching1D& stretch1d,
                                size_t d, level_t& discreteVectorLevel) const;

  /*
   * calculates the left and right neighbor index and levels
   * @param level current level of the node
   * @param index current index of the node
   * @param leftLevel left neighbor level
   * @param leftIndex left neighbor index
   * @param rightLevel right neighbor level
   * @param rightIndex right neighbor index
   */
  void calculateNeighborSpecs(level_t level, index_t index, level_t& leftLevel,
                              index_t& leftIndex, level_t& rightLevel, index_t& rightIndex) const;

  // explicitly hide BoundingBox::unserialize methods to prevent -Woverloaded-virtual warnings
  using BoundingBox::unserialize;
 public:
  /**
   * Constructor for Stretching
   *
   * initializes the Stretching with using no stretching
   *
   * @param dimension number of the dimensions used with the grid
   */
  explicit Stretching(size_t dimension);

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

  /**
   * Destructor.
   */
  ~Stretching() override;

  /*
   * Gets the node position defined in the input parameters
   *
   * @param level level of the node
   * @param index index of the node
   * @param d dimension of the node
   */
  inline double getCoordinate(level_t level, index_t index, size_t d) const {
    if (level == 0) {
      if (index == 0) {
        return boundingBox1Ds[d].leftBoundary;
      } else {
        return boundingBox1Ds[d].rightBoundary;
      }
    } else {
      return stretchingXform(level, index, d);
    }
  }

  /*
   * Returns the type of the Stretching of the dimension given
   *
   * @param d dimension index of the type array
   *
   */
  inline const Stretching1D& getStretching1D(size_t d) const {
    return stretching1Ds[d];
  }

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
  void getAdjacentPositions(level_t level, index_t index, size_t d,
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
  std::vector<level_t> getDiscreteVectorLevel() const;

  void calculateNeighborLookup(level_t maxlevel) const;

  /**
   * Serialize the Stretching into a string.
   *
   * @param version the serialization version of the file
   * @return string that contains all Stretching information
   */
  std::string serialize(int version = SERIALIZATION_VERSION) const override;

  /**
   * Serialize the Stretching into a stream.
   *
   * @param ostream reference to a stream into that all Stretching information is written
   * @param version the serialization version of the file
   */
  void serialize(std::ostream& ostream, int version = SERIALIZATION_VERSION) const override;

  /**
   * Unserialize from a string.
   *
   * @param istr string which contains a serialized Stretching
   * @param mode stretching mode ("analytic" or "discrete")
   * @param version the serialization version of the file
   */
  void unserialize(const std::string& istr, const std::string& mode, int version);

  /**
   * Unserialize from a stream.
   *
   * @param istr stream which contains a serialized Stretching
   * @param mode stretching mode ("analytic" or "discrete")
   * @param version the serialization version of the file
   */
  void unserialize(std::istream& istr, const std::string& mode, int version);
};

}  // namespace base
}  // namespace sgpp

#endif /* STRETCHING_HPP */
