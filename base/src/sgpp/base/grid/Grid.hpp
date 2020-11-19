// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/CoarseningConfiguration.hpp>
#include <sgpp/base/grid/RefinementConfiguration.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/globaldef.hpp>

#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * enum to address different gridtypes in a standardized way
 *
 */
enum class GridType {
  Linear,                              //  0
  LinearStretched,                     //  1
  LinearL0Boundary,                    //  2
  LinearBoundary,                      //  3
  LinearStretchedBoundary,             //  4
  LinearTruncatedBoundary,             //  5
  ModLinear,                           //  6
  Poly,                                //  7
  PolyBoundary,                        //  8
  ModPoly,                             //  9
  ModWavelet,                          // 10
  ModBspline,                          // 11
  Prewavelet,                          // 12
  SquareRoot,                          // 13
  Periodic,                            // 14
  LinearClenshawCurtisBoundary,        // 15
  Bspline,                             // 16
  BsplineBoundary,                     // 17
  BsplineClenshawCurtis,               // 18
  Wavelet,                             // 19
  WaveletBoundary,                     // 20
  FundamentalSpline,                   // 21
  ModFundamentalSpline,                // 22
  ModBsplineClenshawCurtis,            // 23
  LinearStencil,                       // 24
  ModLinearStencil,                    // 25
  PolyClenshawCurtisBoundary,          // 26
  PolyClenshawCurtis,                  // 27
  LinearClenshawCurtis,                // 28
  ModPolyClenshawCurtis,               // 29
  ModLinearClenshawCurtis,             // 30
  NaturalBsplineBoundary,              // 31
  NakBsplineBoundary,                  // 32
  ModNakBspline,                       // 33
  WeaklyFundamentalSplineBoundary,     // 34
  WeaklyFundamentalNakSplineBoundary,  // 35
  ModWeaklyFundamentalNakSpline,       // 36
  FundamentalSplineBoundary,           // 37
  FundamentalNakSplineBoundary,        // 38
};

/**
 * Enum to define all possible grid "super" types (used for GeneralGridConfiguration)
 */
enum class GeneralGridType {
  RegularSparseGrid,
  RefinedCoarsenedSparseGrid,
  ComponentGrid,
  GeometryAwareSparseGrid
};

/**
 * Struct to cluster general grid information for different grid "super" types
 */
struct GeneralGridConfiguration {
  // Grid "super" types
  GeneralGridType generalType_ = GeneralGridType::RegularSparseGrid;
  /// Grid Type, see enum
  sgpp::base::GridType type_ = GridType::Linear;
  /// number of dimensions
  size_t dim_ = 0;
  /// number of levels
  int level_ = 3;
  /// vector of levels for each dimension
  /// TODO(Sebastian Kreisel): initialize with some default value!
  std::vector<size_t> levelVector_;
  /// max. polynomial degree for poly basis
  size_t maxDegree_ = 1;
  /// level of boundary grid
  level_t boundaryLevel_ = 0;
  /// string to serialized grid
  std::string filename_ = "";
  /// subgrid selection value t
  double t_ = 0.0;
  /// virtual destructor, since GeneralGridConfiguration is used as base class
  virtual ~GeneralGridConfiguration() {}
};

/**
 * structure that can be used by applications to cluster regular grid information
 */
struct RegularGridConfiguration : GeneralGridConfiguration {
  RegularGridConfiguration() { generalType_ = GeneralGridType::RegularSparseGrid; }

  ~RegularGridConfiguration() override {}
};

/**
 * Structure that can be used by applications to cluster combi grid information
 */
struct CombiGridConfiguration : GeneralGridConfiguration {
  // The level_ member is replaced by a level vector
  std::vector<int> levels;
  CombiGridConfiguration() {
    generalType_ = GeneralGridType::ComponentGrid;
    level_ = -1;
  }

  ~CombiGridConfiguration() override {}
};

/**
 * Enum that is used to set the type of refinement threshold (percentage/relative based or
 * absolute).
 */
enum class AdaptivityThresholdType {
  Relative,
  Absolute,
};

/**
 * structure that can be used by application to define adaptivity strategies
 */
struct AdaptivityConfiguration {
  /// number of refinements
  size_t numRefinements_ = 1;
  /// threshold type
  AdaptivityThresholdType thresholdType_ = AdaptivityThresholdType::Absolute;
  /// refinement threshold for surpluses
  double refinementThreshold_ = 0.0;
  /// coarsening threshold for surpluses
  double coarseningThreshold_ = 0.0;
  /// prevent coarsening of initial grid points, needed for some decompositions
  bool coarsenInitialPoints_ = false;
  /// refinement type: false: classic, true: maxLevel
  bool maxLevelType_ = false;
  /// max. number of points to be refined
  size_t numRefinementPoints_ = 5;
  /// max. number of points to be coarsened
  size_t numCoarseningPoints_ = 5;
  /// max. percent of points to be refined/coarsened
  double percent_ = 1.0;
  /// other refinement strategy, that is more expensive, but yields better results
  bool errorBasedRefinement_ = false;
  /// threshold for convergence in case error based refinement is applied
  double errorConvergenceThreshold_ = 0.001;
  /// amount of error values to consider when checking for convergence in case of error based
  /// refinement
  size_t errorBufferSize_ = 3;
  /// minimum amount of iterations before the next refinement is allowed to happen in case of error
  /// based refinement
  size_t errorMinInterval_ = 0;
  /// refinement will be triggered each refinementPeriod instances (approximately) in case of non
  /// error based refinement
  size_t refinementPeriod_ = 1;
  /// refinement indicator
  RefinementFunctorType refinementFunctorType_ = RefinementFunctorType::Surplus;
  /// in case of zero corssing based refinement: determines if evaluations should be precomupted
  bool precomputeEvaluations_ = true;
  /// determines if finer grid levels should be penalized when finding points to refine
  bool levelPenalize_ = false;
  /// in case of data based refinements: determines the scaling coefficients for each class
  std::vector<double> scalingCoefficients_ = std::vector<double>();
  /// coarsening indicator
  CoarseningFunctorType coarseningFunctorType_ = CoarseningFunctorType::Surplus;
};

/**
 * abstract base class for all types grids used in sgpp the class gives pure virtual function
 * definitions that have to be implemented by all types of grids
 */
class Grid {
 public:
  /**
   * delete copy constructor
   * @param other
   */
  Grid(const Grid& other) = delete;

  /**
   * creates a grid defined by the grid configuration
   *
   * @param gridConfig grid configuration
   * @return grid
   */
  static Grid* createGrid(RegularGridConfiguration gridConfig);

  /**
   * creates a stencil for a linear grid (without boundaries)
   *
   * <table border="0"><tr>
   * <td>\image html "createLinearGridStencil_C2J-small.png" "Level 4 sparse grid"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createLinearGridStencil(size_t dim);

  /**
   * creates a stencil for a modified linear grid (without boundaries)
   *
   * <table border="0"><tr>
   * <td>\image html "createModLinearGridStencil_C2J-small.png" "Level 4 sparse grid"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createModLinearGridStencil(size_t dim);

  /**
   * Creates and returns a grid without grid points on the boundary (zero boundary conditions) with
   * piecewise linear basis functions
   *
   * <table border="0"><tr>
   * <td>\image html "createLinearGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createLinearGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createLinearGrid(size_t dim);

  /**
   * creates a linear stretched grid without boundaries
   *
   * <table border="0"><tr>
   * <td>\image html "createLinearStretchedGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createLinearStretchedGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createLinearStretchedGrid(size_t dim);

  /**
   * creates a linear boundary grid
   * <table border="0"><tr>
   * <td>\image html "createLinearBoundaryGrid_C2,_0J-small.png"
   * "Level 4 sparse grid with boundaryLevel = 0"</td>
   * <td>\image html "createLinearBoundaryGrid_C2,_1J-small.png"
   * "Level 4 sparse grid with boundaryLevel = 1"</td>
   * <td>\image html "createLinearBoundaryGrid_C2,_2J-small.png"
   * "Level 4 sparse grid with boundaryLevel = 2"</td>
   * <td>\image html "hiba_createLinearGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param boundaryLevel on which level the boundary grid points and
   *                      basis functions should be added;
   *                      the default is 1, which results in a grid with
   *                      the same resolution on the boundary as on the
   *                      main axis
   * @return grid
   */
  static Grid* createLinearBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * creates a linearstretched truncated boundary grid
   *
   * <table border="0"><tr>
   * <td>\image html "createLinearStretchedBoundaryGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createLinearStretchedBoundaryGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   */
  static Grid* createLinearStretchedBoundaryGrid(size_t dim);

  /**
   * creates a linear Clenshaw-Curtis boundary grid
   *
   * @param dim the grid's dimension
   * @param boundaryLevel level of the boundary
   * @return grid
   */
  static Grid* createLinearClenshawCurtisBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * creates a linear Clenshaw-Curtis grid
   *
   * <table border="0"><tr>
   * <td>\image html "createLinearClenshawCurtisGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createLinearClenshawCurtisGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createLinearClenshawCurtisGrid(size_t dim);

  /**
   * creates a modified linear Clenshaw-Curtis grid
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createModLinearClenshawCurtisGrid(size_t dim);

  /**
   * creates a modified linear grid
   *
   * <table border="0"><tr>
   * <td>\image html "createModLinearGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createModLinearGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createModLinearGrid(size_t dim);

  /**
   * creates a polynomial grid
   *
   * <table border="0"><tr>
   * <td>\image html "createPolyGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createPolyGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static Grid* createPolyGrid(size_t dim, size_t degree);

  /**
   * creates a polynomial grid with truncated boundary
   *
   * <table border="0"><tr>
   * <td>\image html "createPolyBoundaryGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createPolyBoundaryGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @param boundaryLevel on which level the boundary grid points and
   *                      basis functions should be added;
   *                      the default is 1, which results in a grid with
   *                      the same resolution on the boundary as on the
   *                      main axis
   * @return grid
   */
  static Grid* createPolyBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel = 1);

  /**
   * creates a poly Clenshaw Curtis boundary grid with clenshaw curtis points
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @param boundaryLevel level at which boundary points are added
   * @return grid
   */
  static Grid* createPolyClenshawCurtisBoundaryGrid(size_t dim, size_t degree,
                                                    level_t boundaryLevel = 1);

  /**
   * creates a poly grid with clenshaw curtis points
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static Grid* createPolyClenshawCurtisGrid(size_t dim, size_t degree);

  /**
   * creates a modified poly grid with clenshaw curtis points
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static Grid* createModPolyClenshawCurtisGrid(size_t dim, size_t degree);

  /**
   * creates a modified polynomial grid
   *
   * <table border="0"><tr>
   * <td>\image html "createModPolyGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createModPolyGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the polynom's max. degree
   * @return grid
   */
  static Grid* createModPolyGrid(size_t dim, size_t degree);

  /**
   * creates a wavelet grid
   *
   * <table border="0"><tr>
   * <td>\image html "createWaveletGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createWaveletGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createWaveletGrid(size_t dim);

  /**
   * creates a wavelet trapezoid boundary grid
   *
   * <table border="0"><tr>
   * <td>\image html "createWaveletBoundaryGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createWaveletBoundaryGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param boundaryLevel on which level the boundary grid points and
   *                      basis functions should be added;
   *                      the default is 1, which results in a grid with
   *                      the same resolution on the boundary as on the
   *                      main axis
   * @return grid
   */
  static Grid* createWaveletBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * creates a modified wavelet grid
   *
   * <table border="0"><tr>
   * <td>\image html "createModWaveletGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createModWaveletGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createModWaveletGrid(size_t dim);

  /**
   * creates a B-spline grid
   *
   * <table border="0"><tr>
   * <td>\image html "createBsplineGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createBsplineGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static Grid* createBsplineGrid(size_t dim, size_t degree);

  /**
   * creates a B-spline trapezoid boundary grid
   *
   * <table border="0"><tr>
   * <td>\image html "createBsplineBoundaryGrid_C2,_3J-small.png" "Level 4* sparse grid"</td>
   * <td>\image html "hiba_createBsplineBoundaryGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @param boundaryLevel on which level the boundary grid points and
   *                      basis functions should be added;
   *                      the default is 1, which results in a grid with
   *                      the same resolution on the boundary as on the
   *                      main axis
   * @return grid
   */
  static Grid* createBsplineBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel = 1);

  /**
   * creates a B-spline Clenshaw-Curtis grid
   *
   * <table border="0"><tr>
   * <td>\image html "createBsplineClenshawCurtisGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createBsplineClenshawCurtisGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @param boundaryLevel on which level the boundary grid points and
   *                      basis functions should be added;
   *                      the default is 1, which results in a grid with
   *                      the same resolution on the boundary as on the
   *                      main axis
   * @return grid
   */
  static Grid* createBsplineClenshawCurtisGrid(size_t dim, size_t degree,
                                               level_t boundaryLevel = 1);

  /**
   * creates a modified B-spline grid
   *
   * <table border="0"><tr>
   * <td>\image html "createModBsplineGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createModBsplineGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static Grid* createModBsplineGrid(size_t dim, size_t degree);

  /**
   * creates a modified B-spline Clenshaw-Curtis grid
   *
   * <table border="0"><tr>
   * <td>\image html "createModBsplineClenshawCurtisGrid_C2,_3J-small.png"
   * "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createModBsplineClenshawCurtisGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static Grid* createModBsplineClenshawCurtisGrid(size_t dim, size_t degree);

  /**
   * creates a fundamental spline grid
   *
   * <table border="0"><tr>
   * <td>\image html "createFundamentalSplineGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createFundamentalSplineGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static Grid* createFundamentalSplineGrid(size_t dim, size_t degree);

  /**
   * creates a modified fundamental spline grid
   *
   * <table border="0"><tr>
   * <td>\image html "createModFundamentalSplineGrid_C2,_3J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createModFundamentalSplineGrid_C2,_3J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @return grid
   */
  static Grid* createModFundamentalSplineGrid(size_t dim, size_t degree);

  /**
   * creates a prewavelet grid
   *
   * <table border="0"><tr>
   * <td>\image html "createPrewaveletGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createPrewaveletGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createPrewaveletGrid(size_t dim);

  /**
   * creates a square root grid (h-grid)
   *
   * <table border="0"><tr>
   * <td>\image html "createSquareRootGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createSquareRootGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createSquareRootGrid(size_t dim);

  /**
   * creates a truncated boundary grid=contains all the gridpoints of the fullgrids which have
   * \f$|l|<level and li>=l_user\f$
   *
   * <table border="0"><tr>
   * <td>\image html "createLinearTruncatedBoundaryGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createLinearTruncatedBoundaryGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createLinearTruncatedBoundaryGrid(size_t dim);

  /**
   * creates a periodic grid
   *
   * <table border="0"><tr>
   * <td>\image html "createPeriodicGrid_C2J-small.png" "Level 4 sparse grid"</td>
   * <td>\image html "hiba_createPeriodicGrid_C2J-small.png"
   * "Hierarchical basis functions up to level 3"</td>
   * </tr></table>
   *
   * @param dim the grid's dimension
   * @return grid
   */
  static Grid* createPeriodicGrid(size_t dim);

  /**
   * creates a not a knot B-Spline boundary grid
   *
   * @param dim the grid's dimension
   * @param degree the B-spline degree
   * @param boundaryLevel the level of the boundary grid
   * @return grid
   */
  static Grid* createNaturalBsplineBoundaryGrid(size_t dim, size_t degree,
                                                level_t boundaryLevel = 1);
  static Grid* createNakBsplineBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel = 1);
  static Grid* createModNakBsplineGrid(size_t dim, size_t degree);
  static Grid* createWeaklyFundamentalSplineBoundaryGrid(size_t dim, size_t degree,
                                                         level_t boundaryLevel = 1);
  static Grid* createWeaklyFundamentalNakSplineBoundaryGrid(size_t dim, size_t degree,
                                                            level_t boundaryLevel = 1);
  static Grid* createModWeaklyFundamentalNakSplineGrid(size_t dim, size_t degree);
  static Grid* createFundamentalSplineBoundaryGrid(size_t dim, size_t degree,
                                                   level_t boundaryLevel = 1);
  static Grid* createFundamentalNakSplineBoundaryGrid(size_t dim, size_t degree,
                                                      level_t boundaryLevel = 1);

  /**
   * reads a grid out of a string
   *
   * @param istr string that contains the grid information
   * @return grid
   */
  static Grid* unserialize(const std::string& istr);

  /**
   * reads a grid out of a stream
   * @param istr inputstream that contains the grid information
   * @return grid
   */
  static Grid* unserialize(std::istream& istr);

 protected:
  /**
   * This constructor creates a new GridStorage out of the stream.
   * For derived classes create an own constructor wich takes a std::istream and calls this
   * function. Add your own static unserialize function and add it in typeMap().
   *
   * @param istr inputstream that contains the grid information
   */
  explicit Grid(std::istream& istr);

  /**
   * Constructor initializing the grid storage with the given dimensionality.
   *
   * @param dim dimensionality
   */
  explicit Grid(size_t dim);

  /**
   * Constructor initializing the grid storage with the given BoundingBox.
   *
   * @param boundingBox BoundingBox of the grid
   */
  explicit Grid(BoundingBox& boundingBox);

  /**
   * Constructor initializing the grid storage with the given Stretching.
   *
   * @param stretching Stretching of the grid
   */
  explicit Grid(Stretching& stretching);

 public:
  /**
   * Desctructor
   */
  virtual ~Grid();

  /**
   * copies a grid
   */
  Grid* clone();

  /**
   * creates an equivalent grid without copying the grid points
   * @param numDims number of dimensions
   */
  Grid* createGridOfEquivalentType(size_t numDims);

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
   * sets the GridStorage's BoundingsBox
   */
  virtual void setBoundingBox(BoundingBox& boundingBox);

  /**
   * sets the GridStorage's Stretching pointer to a Stretching object
   */
  virtual void setStretching(Stretching& stretching);

  /**
   * @return reference to a GridGenerator object
   */
  virtual GridGenerator& getGenerator() = 0;

  /**
   * Returns a string that identifies the grid type uniquely
   *
   * @return string that identifies the grid type uniquely
   */
  virtual sgpp::base::GridType getType() = 0;

  /**
   * Returns a string that identifies the grid type uniquely
   *
   * @return string that identifies the grid type uniquely
   */
  std::string getTypeAsString();

  /**
   * Returns the grid type that corresponds to the actual type but does no boundary treatment
   *
   * @return grid type
   */
  sgpp::base::GridType getZeroBoundaryType();

  /**
   * Returns the Basis class associated with the grid
   *
   * @return Basis class associated with the grid
   */
  virtual SBasis& getBasis() = 0;

  /**
   * Serializes grid to a string.
   * Needed for Python compatibility. Calls serialize(std::ostream&).
   *
   * @param ostr string into which the grid is written
   * @param version the serialization version of the file
   */
  void serialize(std::string& ostr, int version = SERIALIZATION_VERSION);

  /**
   * Serializes the grid.
   * Override if additional information need to be saved.
   * Call base function before writing anything!
   *
   * @param ostr stream to which the grid is written
   * @param version the serialization version of the file
   */
  virtual void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION);

  /**
   * Serializes grid to a string.
   * Needed for Java compatibility.
   *
   * @param version the serialization version of the file
   * @returns string into which the grid is written
   */
  std::string serialize(int version = SERIALIZATION_VERSION);

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
  void insertPoint(size_t dim, unsigned int levels[], unsigned int indices[], bool isLeaf);

  /**
   * Returns the number of dimensions
   * @return the number of dimensions
   */
  size_t getDimension() const;

  /**
   * Returns the number of points on the grid
   * @return the number of points on the grid
   */
  size_t getSize() const;

  /**
   * returns the algorithmic dimensions (the dimensions in which the Up Down operations should be
   * applied)
   *
   * @return the algorithmic dimensions
   */
  std::vector<size_t> getAlgorithmicDimensions();

  /**
   * sets the algorithmic dimensions (the dimensions in which the Up Down operations should be
   * applied)
   *
   * @param newAlgoDims std::vector containing the algorithmic dimensions
   */
  void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

  /**
   * Conversion from string to grid type
   *
   * @param gridType grid type as a string
   * @return actual grid type
   */
  static GridType stringToGridType(const std::string& gridType);

 protected:
  /// GridStorage object of the grid
  GridStorage storage;

  typedef Grid* (*Factory)(std::istream&);
  typedef std::map<std::string, Grid::Factory> factoryMap;
  typedef std::map<sgpp::base::GridType, std::string> gridTypeVerboseMap;

  static Grid* nullFactory(std::istream&);

 private:
  /**
   * This method returns a map with all available grid types for serialization
   *
   * @return a map with all available grid types for serialization
   */
  static factoryMap& typeMap();

  /**
   * This method returns a map with string representation for the type of all available grids.
   * @return a map with string representation for the type of all available grids.
   */
  static gridTypeVerboseMap& typeVerboseMap();
};

}  // namespace base
}  // namespace sgpp
