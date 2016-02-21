// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRID_HPP
#define GRID_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <map>
#include <string>
#include <vector>

namespace SGPP {
namespace base {

/**
 * enum to address different gridtypes in a standardized way
 *
 */
enum class GridType {
  Linear,                       //  0
  LinearStretched,              //  1
  LinearL0Boundary,             //  2
  LinearBoundary,               //  3
  LinearStretchedBoundary,      //  4
  LinearTruncatedBoundary,      //  5
  ModLinear,                    //  6
  Poly,                         //  7
  PolyBoundary,                 //  8
  ModPoly,                      //  9
  ModWavelet,                   // 10
  ModBspline,                   // 11
  Prewavelet,                   // 12
  SquareRoot,                   // 13
  Periodic,                     // 14
  LinearClenshawCurtis,         // 15
  Bspline,                      // 16
  BsplineBoundary,              // 17
  BsplineClenshawCurtis,        // 18
  Wavelet,                      // 19
  WaveletBoundary,              // 20
  FundamentalSpline,            // 21
  ModFundamentalSpline,         // 22
  ModBsplineClenshawCurtis,     // 23
  LinearStencil,                // 24
  ModLinearStencil              // 25
};

/**
 * structure that can be used by applications to cluster regular grid information
 */
struct RegularGridConfiguration {
  /// Grid Type, see enum
  SGPP::base::GridType type_;
  /// number of dimensions
  size_t dim_;
  /// number of levels
  int level_;
  /// max. polynomial degree for poly basis
  size_t maxDegree_;
  /// level of boundary grid
  size_t boundaryLevel_;
};

/**
 * structure that can be used by application to define adaptivity strategies
 */
struct AdpativityConfiguration {
  /// number of refinements
  size_t numRefinements_;
  /// refinement threshold for surpluses
  float_t threshold_;
  /// refinement type: false: classic, true: maxLevel
  bool maxLevelType_;
  /// max. number of points to be refined
  size_t noPoints_;
  /// max. percent of points to be refined
  float_t percent_;
};

/**
 * abstract base class for all types grids used in sgpp
 * the class gives pure virtual function definitions that
 * have to be implemented by all types of grids
 */
class Grid {
 public:
  /**
   * creates a stencil for a linear grid without boundaries
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createLinearGridStencil(size_t dim);

  /**
   * creates a stencil for a modified linear grid (without boundaries)
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createModLinearGridStencil(size_t dim);

  /**
   * creates a linear grid without boundaries
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createLinearGrid(size_t dim);

  /**
   * creates a linear stretched grid without boundaries
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createLinearStretchedGrid(size_t dim);

  /**
   * creates a linear boundary grid
   *
   * @param dim the grid's dimension
   * @param boundaryLevel on which level the boundary grid points and
   *                      basis functions should be added;
   *                      the default is 1, which results in a grid with
   *                      the same resolution on the boundary as on the
   *                      main axis
   * @return grid
   */
  static std::unique_ptr<Grid> createLinearBoundaryGrid(size_t dim,
                                        level_t boundaryLevel = 1);

  /**
   * creates a linearstretched truncated boundary grid
   *
   * @param dim the grid's dimension
   */
  static std::unique_ptr<Grid> createLinearStretchedBoundaryGrid(size_t dim);

  /**
   * creates a linear Clenshaw-Curtis grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createLinearClenshawCurtisGrid(size_t dim);

  /**
   * creates a mod linear grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createModLinearGrid(size_t dim);

  /**
   * creates a polynomial grid
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static std::unique_ptr<Grid> createPolyGrid(size_t dim, size_t degree);

  /**
   * creates a polynomial grid with truncated boundary
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static std::unique_ptr<Grid> createPolyBoundaryGrid(size_t dim, size_t degree);

  /**
   * creates a poly grid
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static std::unique_ptr<Grid> createModPolyGrid(size_t dim, size_t degree);

  /**
   * creates a wavelet grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createWaveletGrid(size_t dim);

  /**
   * creates a wavelet trapezoid boundary grid
   *
   * @param dim the grid's dimension
   */
  static std::unique_ptr<Grid> createWaveletBoundaryGrid(size_t dim);

  /**
   * creates a mod wavelet grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createModWaveletGrid(size_t dim);

  /**
   * creates a Bspline grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createBsplineGrid(size_t dim, size_t degree);

  /**
   * creates a Bspline trapezoid boundary grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createBsplineBoundaryGrid(size_t dim, size_t degree);

  /**
   * creates a Bspline Clenshaw-Curtis grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createBsplineClenshawCurtisGrid(size_t dim, size_t degree);

  /**
   * creates a mod-Bspline grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createModBsplineGrid(size_t dim, size_t degree);

  /**
   * creates a mod-Bspline Clenshaw-Curtis grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createModBsplineClenshawCurtisGrid(size_t dim, size_t degree);

  /**
   * creates a fundamental spline grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createFundamentalSplineGrid(size_t dim, size_t degree);

  /**
   * creates a mod-fundamental spline grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static std::unique_ptr<Grid> createModFundamentalSplineGrid(size_t dim, size_t degree);

  /**
   * creates a prewavelet grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createPrewaveletGrid(size_t dim);

  /**
   * creates a square root grid(h-grid)
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createSquareRootGrid(size_t dim);

  /**
   * creates a truncated boundary grid=contains all the gridpoints of the fullgrids which have \f$|l|<level and li>=l_user\f$
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createLinearTruncatedBoundaryGrid(size_t dim);

  /**
   * creates a periodic grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static std::unique_ptr<Grid> createPeriodicGrid(size_t dim);

  /**
   * reads a grid out of a string
   *
   * @param istr string that contains the grid information
   * @return grid
   */
  static std::unique_ptr<Grid> unserialize(const std::string& istr);

  /**
   * reads a grid out of a stream
   * @param istr inputstream that contains the grid information
   * @return grid
   */
  static std::unique_ptr<Grid> unserialize(std::istream& istr);

 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   * For derived classes create an own constructor wich takes a std::istream and calls
   * this function. Add your own static unserialize function and add it in typeMap().
   *
   * @param istr inputstream that contains the grid information
   */
  explicit Grid(std::istream& istr);

  /**
   * Constructor initializing the grid storage with the given
   * dimensionality.
   *
   * @param dim dimensionality
   */
  explicit Grid(size_t dim);

  /**
   * Constructor initializing the grid storage with the given
   * BoundingBox.
   *
   * @param BB BoundingBox of the grid
   */
  explicit Grid(BoundingBox& BB);

  /**
   * Constructor initializing the grid storage with the given
   * Stretching.
   *
   * @param BB Stretching of the grid
   */
  explicit Grid(Stretching& BB);

 public:
  /**
   * Desctructor
   */
  virtual ~Grid();

  /**
   * gets a reference to the GridStorage object
   *
   * @return reference to the GridStorage obeject
   */
  virtual GridStorage& getStorage();

  /**
   * gets a reference to the GridStorage's BoundingsBox object
   *
   * @return reference to the GridStorage's BoundingsBox object
   */
  virtual BoundingBox& getBoundingBox();

  /**
   * gets a reference to the GridStorage's Stretching object
   *
   * @return reference to the GridStorage's Stretching object
   */
  virtual Stretching& getStretching();

  /**
   * sets the GridStorage's BoundingsBox pointer to a BoundingBox object
   *
   * @return pointer to the GridStorage's BoundingsBox object
   */
  virtual void setBoundingBox(BoundingBox& bb);

  /**
   * sets the GridStorage's Stretching pointer to a Stretching object
   *
   * @return pointer to the GridStorage's Stretching object
   */
  virtual void setStretching(Stretching& bb);

  /**
   * @return reference to a GridGenerator object
   */
  virtual GridGenerator& getGenerator() = 0;

  /**
   * Returns a string that identifies the grid type uniquely
   *
   * @return string that identifies the grid type uniquely
   */
  virtual SGPP::base::GridType getType() = 0;

  /**
   * Returns the Basis class associated with the grid
   *
   * @return Basis class associated with the grid
   */
  virtual const SBasis& getBasis() = 0;

  /**
   * Serializes grid to a string.
   * Needed for Python compatibility. Calls serialize(std::ostream&).
   *
   * @param ostr string into which the grid is written
   */
  void serialize(std::string& ostr);

  /**
   * Serializes the grid.
   * Override if additional information need to be saved.
   * Call base function before writing anything!
   *
   * @param ostr stream to which the grid is written
   */
  virtual void serialize(std::ostream& ostr);

  /**
   * Serializes grid to a string.
   * Needed for Java compatibility.
   *
   * @returns string into which the grid is written
   */
  std::string serialize();

  /**
   * Refine grid
   * Refine the given number of points on the grid according to the vector
   *
   * @param vector DataVector vector with errors for each basis function or alpha-vector
   * @param numOfPoints integer number of points to refine
   */
  void refine(DataVector& vector, int numOfPoints);

  /**
   * Insert one point to the grid
   *
   * @param dim dimension of the grid
   * @param levels array with levels of the point
   * @param indices array with indices of the point
   * @param isLeaf indicator whether the point is a leaf
   */
  void insertPoint(size_t dim, unsigned int levels[], unsigned int indices[],
                   bool isLeaf);

  /**
   * Returns the number of points on the grid
   * @return the number of points on the grid
   */
  size_t getSize();

  /**
   * returns the algorithmic dimensions (the dimensions in which the Up Down
   * operations should be applied)
   *
   * @return the algorithmic dimensions
   */
  std::vector<size_t> getAlgorithmicDimensions();

  /**
   * sets the algorithmic dimensions (the dimensions in which the Up Down
   * operations should be applied)
   *
   * @param newAlgoDims std::vector containing the algorithmic dimensions
   */
  void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

 protected:
  /// GridStorage object of the grid
  GridStorage storage;

  typedef std::unique_ptr<Grid> (*Factory)(std::istream&);
  typedef std::map<std::string, Grid::Factory> factoryMap;
  typedef std::map<SGPP::base::GridType, std::string> gridTypeVerboseMap;

  static std::unique_ptr<Grid> nullFactory(std::istream&);

 private:
  /**
   * This method returns a map with all available grid types for serialization
   *
   * @return a map with all available grid types for serialization
   */
  static factoryMap& typeMap();

  /**
   *  This method returns a map with string representation for the type
   *  of all available grids.
   * @return a map with string representation for the type
   *  of all available grids.
   */
  static gridTypeVerboseMap& typeVerboseMap();
};

}  // namespace base
}  // namespace SGPP

#endif /* GRID_HPP */
