// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/tools/sle/system/CloneableSLE.hpp>

#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>

#include <cstddef>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace sgpp {
namespace base {

/**
 * Linear system of the hierarchization in a sparse grid.
 */
class HierarchisationSLE : public CloneableSLE {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param grid              sparse grid
   */
  explicit HierarchisationSLE(Grid& grid) : HierarchisationSLE(grid, grid.getStorage()) {}

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param grid              sparse grid
   * @param gridStorage       custom grid storage (use basis function
   *                          according to grid, but use another set of
   *                          grid points according to gridStorage)
   */
  HierarchisationSLE(Grid& grid, GridStorage& gridStorage)
      : CloneableSLE(), grid(grid), gridStorage(gridStorage), basisType(INVALID) {
    // initialize the correct basis (according to the grid)
    if (grid.getType() == GridType::Bspline) {
      bsplineBasis = std::unique_ptr<SBsplineBase>(
          new SBsplineBase(dynamic_cast<BsplineGrid&>(grid).getDegree()));
      basisType = BSPLINE;
    } else if (grid.getType() == GridType::BsplineBoundary) {
      bsplineBoundaryBasis = std::unique_ptr<SBsplineBoundaryBase>(
          new SBsplineBoundaryBase(dynamic_cast<BsplineBoundaryGrid&>(grid).getDegree()));
      basisType = BSPLINE_BOUNDARY;
    } else if (grid.getType() == GridType::BsplineClenshawCurtis) {
      bsplineClenshawCurtisBasis =
          std::unique_ptr<SBsplineClenshawCurtisBase>(new SBsplineClenshawCurtisBase(
              dynamic_cast<BsplineClenshawCurtisGrid&>(grid).getDegree()));
      basisType = BSPLINE_CLENSHAW_CURTIS;
    } else if (grid.getType() == GridType::ModBspline) {
      modBsplineBasis = std::unique_ptr<SBsplineModifiedBase>(
          new SBsplineModifiedBase(dynamic_cast<ModBsplineGrid&>(grid).getDegree()));
      basisType = BSPLINE_MODIFIED;
    } else if (grid.getType() == GridType::ModBsplineClenshawCurtis) {
      modBsplineClenshawCurtisBasis = std::unique_ptr<SBsplineModifiedClenshawCurtisBase>(
          new SBsplineModifiedClenshawCurtisBase(
              dynamic_cast<ModBsplineClenshawCurtisGrid&>(grid).getDegree()));
      basisType = BSPLINE_MODIFIED_CLENSHAW_CURTIS;
    } else if (grid.getType() == GridType::FundamentalSpline) {
      fundamentalSplineBasis = std::unique_ptr<SFundamentalSplineBase>(
          new SFundamentalSplineBase(dynamic_cast<FundamentalSplineGrid&>(grid).getDegree()));
      basisType = FUNDAMENTAL_SPLINE;
    } else if (grid.getType() == GridType::ModFundamentalSpline) {
      modFundamentalSplineBasis =
          std::unique_ptr<SFundamentalSplineModifiedBase>(new SFundamentalSplineModifiedBase(
              dynamic_cast<ModFundamentalSplineGrid&>(grid).getDegree()));
      basisType = FUNDAMENTAL_SPLINE_MODIFIED;
    } else if (grid.getType() == GridType::Linear) {
      linearBasis = std::unique_ptr<SLinearBase>(new SLinearBase());
      basisType = LINEAR;
    } else if (grid.getType() == GridType::LinearBoundary) {
      linearL0BoundaryBasis = std::unique_ptr<SLinearBoundaryBase>(new SLinearBoundaryBase());
      basisType = LINEAR_BOUNDARY;
    } else if (grid.getType() == GridType::LinearClenshawCurtis) {
      linearClenshawCurtisBasis =
          std::unique_ptr<SLinearClenshawCurtisBase>(new SLinearClenshawCurtisBase());
      basisType = LINEAR_CLENSHAW_CURTIS;
    } else if (grid.getType() == GridType::LinearClenshawCurtisBoundary) {
      linearClenshawCurtisBoundaryBasis = std::unique_ptr<SLinearClenshawCurtisBoundaryBase>(
          new SLinearClenshawCurtisBoundaryBase());
      basisType = LINEAR_CLENSHAW_CURTIS_BOUNDARY;
    } else if (grid.getType() == GridType::ModLinear) {
      modLinearBasis = std::unique_ptr<SLinearModifiedBase>(new SLinearModifiedBase());
      basisType = LINEAR_MODIFIED;
    } else if (grid.getType() == GridType::Wavelet) {
      waveletBasis = std::unique_ptr<SWaveletBase>(new SWaveletBase());
      basisType = WAVELET;
    } else if (grid.getType() == GridType::WaveletBoundary) {
      waveletBoundaryBasis = std::unique_ptr<SWaveletBoundaryBase>(new SWaveletBoundaryBase());
      basisType = WAVELET_BOUNDARY;
    } else if (grid.getType() == GridType::ModWavelet) {
      modWaveletBasis = std::unique_ptr<SWaveletModifiedBase>(new SWaveletModifiedBase());
      basisType = WAVELET_MODIFIED;
    } else if (grid.getType() == GridType::NakBsplineModified) {
      nakBsplineModifiedBasis =
          std::unique_ptr<base::SNakBsplineModifiedBase>(new SNakBsplineModifiedBase(
              dynamic_cast<base::NakBsplineModifiedGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINE_MODIFIED;
    } else if (grid.getType() == GridType::ModPoly) {
      modPolyBasis = std::unique_ptr<SPolyModifiedBase>(
          new SPolyModifiedBase(dynamic_cast<ModPolyGrid&>(grid).getDegree()));
      basisType = MOD_POLY;
    } else if (grid.getType() == GridType::Poly) {
      polyBasis =
          std::unique_ptr<SPolyBase>(new SPolyBase(dynamic_cast<PolyGrid&>(grid).getDegree()));
      basisType = POLY;
    } else if (grid.getType() == GridType::PolyBoundary) {
      polyBoundaryBasis = std::unique_ptr<SPolyBoundaryBase>(
          new SPolyBoundaryBase(dynamic_cast<base::PolyBoundaryGrid&>(grid).getDegree()));
      basisType = POLYBOUNDARY;
    } else if (grid.getType() == GridType::NakBspline) {
      nakBsplineBasis = std::unique_ptr<SNakBsplineBase>(
          new SNakBsplineBase(dynamic_cast<NakBsplineGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINE;
    } else if (grid.getType() == GridType::NakBsplineBoundary) {
      nakBsplineBoundaryBasis = std::unique_ptr<SNakBsplineBoundaryBase>(
          new SNakBsplineBoundaryBase(dynamic_cast<NakBsplineBoundaryGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINEBOUNDARY;
    } else if (grid.getType() == GridType::NakBsplineExtended) {
      nakBsplineExtendedBasis =
          std::unique_ptr<SNakBsplineExtendedBase>(new base::SNakBsplineExtendedBase(
              dynamic_cast<NakBsplineExtendedGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINE_EXTENDED;
    } else {
      throw std::invalid_argument("HierarchisationSLE: Grid type not supported.");
    }
  }

  /**
   * @param i     row index
   * @param j     column index
   * @return      whether the i-th grid point lies in the support of
   *              the j-th basis function
   */
  inline bool isMatrixEntryNonZero(size_t i, size_t j) override {
    return (evalBasisFunctionAtGridPoint(j, i) != 0.0);
  }

  /**
   * @param i     row index
   * @param j     column index
   * @return      value of the j-th basis function at the i-th grid point
   */
  inline double getMatrixEntry(size_t i, size_t j) override {
    return evalBasisFunctionAtGridPoint(j, i);
  }

  /**
   * @return          sparse grid
   */
  Grid& getGrid() { return grid; }

  /**
   * @return              grid storage
   */
  GridStorage& getGridStorage() { return gridStorage; }

  size_t getDimension() const override { return gridStorage.getSize(); }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<CloneableSLE>& clone) const override {
    clone = std::unique_ptr<CloneableSLE>(new HierarchisationSLE(grid, gridStorage));
  }

 protected:
  /// sparse grid
  Grid& grid;
  /// grid storage
  GridStorage& gridStorage;

  /// B-spline basis
  std::unique_ptr<SBsplineBase> bsplineBasis;
  /// B-spline boundary basis
  std::unique_ptr<SBsplineBoundaryBase> bsplineBoundaryBasis;
  /// B-spline Clenshaw-Curtis basis
  std::unique_ptr<SBsplineClenshawCurtisBase> bsplineClenshawCurtisBasis;
  /// modified B-spline basis
  std::unique_ptr<SBsplineModifiedBase> modBsplineBasis;
  /// modified B-spline Clenshaw-Curtis basis
  std::unique_ptr<SBsplineModifiedClenshawCurtisBase> modBsplineClenshawCurtisBasis;
  /// fundamental spline basis
  std::unique_ptr<SFundamentalSplineBase> fundamentalSplineBasis;
  /// modified fundamental spline basis
  std::unique_ptr<SFundamentalSplineModifiedBase> modFundamentalSplineBasis;
  /// linear basis
  std::unique_ptr<SLinearBase> linearBasis;
  /// linear boundary basis
  std::unique_ptr<SLinearBoundaryBase> linearL0BoundaryBasis;
  /// linear Clenshaw-Curtis basis
  std::unique_ptr<SLinearClenshawCurtisBase> linearClenshawCurtisBasis;
  /// linear Clenshaw-Curtis boundary basis
  std::unique_ptr<SLinearClenshawCurtisBoundaryBase> linearClenshawCurtisBoundaryBasis;
  /// modified linear basis
  std::unique_ptr<SLinearModifiedBase> modLinearBasis;
  /// wavelet basis
  std::unique_ptr<SWaveletBase> waveletBasis;
  /// wavelet boundary basis
  std::unique_ptr<SWaveletBoundaryBase> waveletBoundaryBasis;
  /// modified wavelet basis
  std::unique_ptr<SWaveletModifiedBase> modWaveletBasis;
  /// not-a-knot B-spline Boundary basis
  std::unique_ptr<SNakBsplineBoundaryBase> nakBsplineBoundaryBasis;
  /// not-a-knot B-spline modified basis
  std::unique_ptr<SNakBsplineModifiedBase> nakBsplineModifiedBasis;
  /// mod poly basis
  std::unique_ptr<SPolyModifiedBase> modPolyBasis;
  /// poly basis
  std::unique_ptr<SPolyBase> polyBasis;
  /// poly boundary basis
  std::unique_ptr<SPolyBoundaryBase> polyBoundaryBasis;
  /// nak Bspline basis
  std::unique_ptr<SNakBsplineBase> nakBsplineBasis;
  /// nak Bspline extended basis
  std::unique_ptr<SNakBsplineExtendedBase> nakBsplineExtendedBasis;

  /// type of grid/basis functions
  enum {
    INVALID,
    BSPLINE,
    BSPLINE_BOUNDARY,
    BSPLINE_CLENSHAW_CURTIS,
    BSPLINE_MODIFIED,
    BSPLINE_MODIFIED_CLENSHAW_CURTIS,
    FUNDAMENTAL_SPLINE,
    FUNDAMENTAL_SPLINE_MODIFIED,
    LINEAR,
    LINEAR_BOUNDARY,
    LINEAR_CLENSHAW_CURTIS,
    LINEAR_CLENSHAW_CURTIS_BOUNDARY,
    LINEAR_MODIFIED,
    WAVELET,
    WAVELET_BOUNDARY,
    WAVELET_MODIFIED,
    NAK_BSPLINEBOUNDARY,
    NAK_BSPLINE_MODIFIED,
    MOD_POLY,
    POLY,
    POLYBOUNDARY,
    NAK_BSPLINE,
    NAK_BSPLINE_EXTENDED
  } basisType;

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th basis function at the
   *                  pointJ-th grid point
   */
  inline double evalBasisFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    if (basisType == BSPLINE) {
      return evalBsplineFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == BSPLINE_BOUNDARY) {
      return evalBsplineBoundaryFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == BSPLINE_CLENSHAW_CURTIS) {
      return evalBsplineClenshawCurtisFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == BSPLINE_MODIFIED) {
      return evalBsplineModifiedFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == BSPLINE_MODIFIED_CLENSHAW_CURTIS) {
      return evalBsplineModifiedClenshawCurtisFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == FUNDAMENTAL_SPLINE) {
      return evalFundamentalSplineFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == FUNDAMENTAL_SPLINE_MODIFIED) {
      return evalFundamentalSplineModifiedFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == LINEAR) {
      return evalLinearFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == LINEAR_BOUNDARY) {
      return evalLinearBoundaryFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == LINEAR_CLENSHAW_CURTIS) {
      return evalLinearClenshawCurtisFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == LINEAR_CLENSHAW_CURTIS_BOUNDARY) {
      return evalLinearClenshawCurtisBoundaryFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == LINEAR_MODIFIED) {
      return evalLinearModifiedFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == WAVELET) {
      return evalWaveletFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == WAVELET_BOUNDARY) {
      return evalWaveletBoundaryFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == WAVELET_MODIFIED) {
      return evalWaveletModifiedFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == NAK_BSPLINEBOUNDARY) {
      return evalNakBsplineBoundaryFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == NAK_BSPLINE_MODIFIED) {
      return evalNakBsplineModifiedFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == MOD_POLY) {
      return evalModPolyFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == POLY) {
      return evalPolyFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == POLYBOUNDARY) {
      return evalPolyBoundaryFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == NAK_BSPLINE) {
      return evalNakBsplineFunctionAtGridPoint(basisI, pointJ);
    } else if (basisType == NAK_BSPLINE_EXTENDED) {
      return evalNakBsplineExtendedFunctionAtGridPoint(basisI, pointJ);
    } else {
      return 0.0;
    }
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th B-spline basis function
   *                  at the pointJ-th grid point
   */
  inline double evalBsplineFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = bsplineBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                 gridStorage.getPointIndex(basisI, t),
                                                 gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th B-spline boundary
   *                  basis function at the pointJ-th grid point
   */
  inline double evalBsplineBoundaryFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = bsplineBoundaryBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                         gridStorage.getPointIndex(basisI, t),
                                                         gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th B-spline Clenshaw-Curtis
   *                  basis function at the pointJ-th grid point
   */
  inline double evalBsplineClenshawCurtisFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = bsplineClenshawCurtisBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getUnitPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th modified B-spline
   *                  basis function at the pointJ-th grid point
   */
  inline double evalBsplineModifiedFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = modBsplineBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                    gridStorage.getPointIndex(basisI, t),
                                                    gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th modified Clenshaw-Curtis
   *                  B-spline basis function at the pointJ-th grid point
   */
  inline double evalBsplineModifiedClenshawCurtisFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = modBsplineClenshawCurtisBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getUnitPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th fundamental spline basis
   *                  function at the pointJ-th grid point
   */
  inline double evalFundamentalSplineFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      if (gridStorage.getPointLevel(pointJ, t) < gridStorage.getPointLevel(basisI, t)) {
        return 0.0;
      } else if (gridStorage.getPointLevel(pointJ, t) == gridStorage.getPointLevel(basisI, t)) {
        if (gridStorage.getPointIndex(pointJ, t) != gridStorage.getPointIndex(basisI, t)) {
          return 0.0;
        }
      } else {
        const double result1d = fundamentalSplineBasis->eval(
            gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
            gridStorage.getPointCoordinate(pointJ, t));

        if (result1d == 0.0) {
          return 0.0;
        }

        result *= result1d;
      }
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th modified fundamental spline
   *                  basis function at the pointJ-th grid point
   */
  inline double evalFundamentalSplineModifiedFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      if (gridStorage.getPointLevel(pointJ, t) < gridStorage.getPointLevel(basisI, t)) {
        return 0.0;
      } else if (gridStorage.getPointLevel(pointJ, t) == gridStorage.getPointLevel(basisI, t)) {
        if (gridStorage.getPointIndex(pointJ, t) != gridStorage.getPointIndex(basisI, t)) {
          return 0.0;
        }
      } else {
        const double result1d = modFundamentalSplineBasis->eval(
            gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
            gridStorage.getPointCoordinate(pointJ, t));

        if (result1d == 0.0) {
          return 0.0;
        }

        result *= result1d;
      }
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th linear
   *                  basis function at the pointJ-th grid point
   */
  inline double evalLinearFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = linearBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                gridStorage.getPointIndex(basisI, t),
                                                gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th linear boundary
   *                  basis function at the pointJ-th grid point
   */
  inline double evalLinearBoundaryFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = linearL0BoundaryBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th linear Clenshaw-Curtis
   *                  basis function at the pointJ-th grid point
   */
  inline double evalLinearClenshawCurtisFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = linearClenshawCurtisBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getUnitPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th linear Clenshaw-Curtis
   *                  boundary basis function at the pointJ-th grid point
   */
  inline double evalLinearClenshawCurtisBoundaryFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = linearClenshawCurtisBoundaryBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getUnitPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th modified linear
   *                  basis function at the pointJ-th grid point
   */
  inline double evalLinearModifiedFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = modLinearBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                   gridStorage.getPointIndex(basisI, t),
                                                   gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th wavelet
   *                  basis function at the pointJ-th grid point
   */
  inline double evalWaveletFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = waveletBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                 gridStorage.getPointIndex(basisI, t),
                                                 gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th wavelet boundary
   *                  basis function at the pointJ-th grid point
   */
  inline double evalWaveletBoundaryFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = waveletBoundaryBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                         gridStorage.getPointIndex(basisI, t),
                                                         gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th modified wavelet
   *                  basis function at the pointJ-th grid point
   */
  inline double evalWaveletModifiedFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = modWaveletBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                    gridStorage.getPointIndex(basisI, t),
                                                    gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th not-a-knot B-spline basis function
   *                  at the pointJ-th grid point
   */
  inline double evalNakBsplineBoundaryFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = nakBsplineBoundaryBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th not-a-knot B-spline modified basis function
   *                  at the pointJ-th grid point
   */
  inline double evalNakBsplineModifiedFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = nakBsplineModifiedBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th mod poly (modified Bungartz polynomials) basis function
   *                  at the pointJ-th grid point
   */
  inline double evalModPolyFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = modPolyBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                 gridStorage.getPointIndex(basisI, t),
                                                 gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th poly basis (Bungartz polynomials) function
   *                  at the pointJ-th grid point
   */
  inline double evalPolyFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = polyBasis->eval(gridStorage.getPointLevel(basisI, t),
                                              gridStorage.getPointIndex(basisI, t),
                                              gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th poly boundary basis (Bungartz polynomialswith with
   * boundary) function at the pointJ-th grid point
   */
  inline double evalPolyBoundaryFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = polyBoundaryBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                      gridStorage.getPointIndex(basisI, t),
                                                      gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th nak Bspline basis function
   *                  at the pointJ-th grid point
   */
  inline double evalNakBsplineFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = nakBsplineBasis->eval(gridStorage.getPointLevel(basisI, t),
                                                    gridStorage.getPointIndex(basisI, t),
                                                    gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }
  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th nak Bspline extended basis function
   *                  at the pointJ-th grid point
   */
  inline double evalNakBsplineExtendedFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = nakBsplineExtendedBasis->eval(
          gridStorage.getPointLevel(basisI, t), gridStorage.getPointIndex(basisI, t),
          gridStorage.getPointCoordinate(pointJ, t));

      if (result1d == 0.0) {
        return 0.0;
      }

      result *= result1d;
    }

    return result;
  }
};
}  // namespace base
}  // namespace sgpp
