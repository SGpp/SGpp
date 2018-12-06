// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATIONSLE_HPP
#define SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATIONSLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>

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
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryCombigridBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>

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
#include <sgpp/base/grid/type/NakBsplineBoundaryCombigridGrid.hpp>
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
namespace optimization {

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
  explicit HierarchisationSLE(base::Grid& grid) : HierarchisationSLE(grid, grid.getStorage()) {
    if (gridStorage.hasBoundingBoxOrStretching()) {
      std::cerr
          << "HierarchisationSLE: Methods do not take bounding box or stretching into account!\n";
    }
  }

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param grid              sparse grid
   * @param gridStorage       custom grid storage (use basis function
   *                          according to grid, but use another set of
   *                          grid points according to gridStorage)
   */
  HierarchisationSLE(base::Grid& grid, base::GridStorage& gridStorage)
      : CloneableSLE(), grid(grid), gridStorage(gridStorage), basisType(INVALID) {
    if (gridStorage.hasBoundingBoxOrStretching()) {
      std::cerr
          << "HierarchisationSLE: Methods do not take bounding box or stretching into account!\n";
    }
    // initialize the correct basis (according to the grid)
    if (grid.getType() == base::GridType::Bspline) {
      bsplineBasis = std::unique_ptr<base::SBsplineBase>(
          new base::SBsplineBase(dynamic_cast<base::BsplineGrid&>(grid).getDegree()));
      basisType = BSPLINE;
    } else if (grid.getType() == base::GridType::BsplineBoundary) {
      bsplineBoundaryBasis =
          std::unique_ptr<base::SBsplineBoundaryBase>(new base::SBsplineBoundaryBase(
              dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree()));
      basisType = BSPLINE_BOUNDARY;
    } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
      bsplineClenshawCurtisBasis =
          std::unique_ptr<base::SBsplineClenshawCurtisBase>(new base::SBsplineClenshawCurtisBase(
              dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree()));
      basisType = BSPLINE_CLENSHAW_CURTIS;
    } else if (grid.getType() == base::GridType::ModBspline) {
      modBsplineBasis = std::unique_ptr<base::SBsplineModifiedBase>(
          new base::SBsplineModifiedBase(dynamic_cast<base::ModBsplineGrid&>(grid).getDegree()));
      basisType = BSPLINE_MODIFIED;
    } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
      modBsplineClenshawCurtisBasis = std::unique_ptr<base::SBsplineModifiedClenshawCurtisBase>(
          new base::SBsplineModifiedClenshawCurtisBase(
              dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid).getDegree()));
      basisType = BSPLINE_MODIFIED_CLENSHAW_CURTIS;
    } else if (grid.getType() == base::GridType::FundamentalSpline) {
      fundamentalSplineBasis =
          std::unique_ptr<base::SFundamentalSplineBase>(new base::SFundamentalSplineBase(
              dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree()));
      basisType = FUNDAMENTAL_SPLINE;
    } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
      modFundamentalSplineBasis = std::unique_ptr<base::SFundamentalSplineModifiedBase>(
          new base::SFundamentalSplineModifiedBase(
              dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree()));
      basisType = FUNDAMENTAL_SPLINE_MODIFIED;
    } else if (grid.getType() == base::GridType::Linear) {
      linearBasis = std::unique_ptr<base::SLinearBase>(new base::SLinearBase());
      basisType = LINEAR;
    } else if (grid.getType() == base::GridType::LinearBoundary) {
      linearL0BoundaryBasis =
          std::unique_ptr<base::SLinearBoundaryBase>(new base::SLinearBoundaryBase());
      basisType = LINEAR_BOUNDARY;
    } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
      linearClenshawCurtisBasis =
          std::unique_ptr<base::SLinearClenshawCurtisBase>(new base::SLinearClenshawCurtisBase());
      basisType = LINEAR_CLENSHAW_CURTIS;
    } else if (grid.getType() == base::GridType::LinearClenshawCurtisBoundary) {
      linearClenshawCurtisBoundaryBasis = std::unique_ptr<base::SLinearClenshawCurtisBoundaryBase>(
          new base::SLinearClenshawCurtisBoundaryBase());
      basisType = LINEAR_CLENSHAW_CURTIS_BOUNDARY;
    } else if (grid.getType() == base::GridType::ModLinear) {
      modLinearBasis = std::unique_ptr<base::SLinearModifiedBase>(new base::SLinearModifiedBase());
      basisType = LINEAR_MODIFIED;
    } else if (grid.getType() == base::GridType::Wavelet) {
      waveletBasis = std::unique_ptr<base::SWaveletBase>(new base::SWaveletBase());
      basisType = WAVELET;
    } else if (grid.getType() == base::GridType::WaveletBoundary) {
      waveletBoundaryBasis =
          std::unique_ptr<base::SWaveletBoundaryBase>(new base::SWaveletBoundaryBase());
      basisType = WAVELET_BOUNDARY;
    } else if (grid.getType() == base::GridType::ModWavelet) {
      modWaveletBasis =
          std::unique_ptr<base::SWaveletModifiedBase>(new base::SWaveletModifiedBase());
      basisType = WAVELET_MODIFIED;
    } else if (grid.getType() == base::GridType::NakBsplineBoundary) {
      nakBsplineBoundaryBasis =
          std::unique_ptr<base::SNakBsplineBoundaryBase>(new base::SNakBsplineBoundaryBase(
              dynamic_cast<base::NakBsplineBoundaryGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINEBOUNDARY;
    } else if (grid.getType() == base::GridType::NakBsplineBoundaryCombigrid) {
      nakBsplineBoundaryCombigridBasis = std::unique_ptr<base::SNakBsplineBoundaryCombigridBase>(
          new base::SNakBsplineBoundaryCombigridBase(
              dynamic_cast<base::NakBsplineBoundaryCombigridGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINEBOUNDARY_COMBIGRID;
    } else if (grid.getType() == base::GridType::NakBsplineModified) {
      nakBsplineModifiedBasis =
          std::unique_ptr<base::SNakBsplineModifiedBase>(new base::SNakBsplineModifiedBase(
              dynamic_cast<base::NakBsplineModifiedGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINE_MODIFIED;
    } else if (grid.getType() == base::GridType::ModPoly) {
      modPolyBasis = std::unique_ptr<base::SPolyModifiedBase>(
          new base::SPolyModifiedBase(dynamic_cast<base::ModPolyGrid&>(grid).getDegree()));
      basisType = MOD_POLY;
    } else if (grid.getType() == base::GridType::Poly) {
      polyBasis = std::unique_ptr<base::SPolyBase>(
          new base::SPolyBase(dynamic_cast<base::PolyGrid&>(grid).getDegree()));
      basisType = POLY;
    } else if (grid.getType() == base::GridType::PolyBoundary) {
      polyBoundaryBasis = std::unique_ptr<base::SPolyBoundaryBase>(
          new base::SPolyBoundaryBase(dynamic_cast<base::PolyBoundaryGrid&>(grid).getDegree()));
      basisType = POLYBOUNDARY;
    } else if (grid.getType() == base::GridType::NakBspline) {
      nakBsplineBasis = std::unique_ptr<base::SNakBsplineBase>(
          new base::SNakBsplineBase(dynamic_cast<base::NakBsplineGrid&>(grid).getDegree()));
      basisType = NAK_BSPLINE;
    } else if (grid.getType() == base::GridType::NakBsplineExtended) {
      nakBsplineExtendedBasis =
          std::unique_ptr<base::SNakBsplineExtendedBase>(new base::SNakBsplineExtendedBase(
              dynamic_cast<base::NakBsplineExtendedGrid&>(grid).getDegree()));
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
  base::Grid& getGrid() { return grid; }

  /**
   * @return              grid storage
   */
  base::GridStorage& getGridStorage() { return gridStorage; }

  size_t getDimension() const override { return gridStorage.getSize(); }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<CloneableSLE>& clone) const override {
    clone = std::unique_ptr<CloneableSLE>(new HierarchisationSLE(grid, gridStorage));
  }

 protected:
  /// sparse grid
  base::Grid& grid;
  /// grid storage
  base::GridStorage& gridStorage;

  /// B-spline basis
  std::unique_ptr<base::SBsplineBase> bsplineBasis;
  /// B-spline boundary basis
  std::unique_ptr<base::SBsplineBoundaryBase> bsplineBoundaryBasis;
  /// B-spline Clenshaw-Curtis basis
  std::unique_ptr<base::SBsplineClenshawCurtisBase> bsplineClenshawCurtisBasis;
  /// modified B-spline basis
  std::unique_ptr<base::SBsplineModifiedBase> modBsplineBasis;
  /// modified B-spline Clenshaw-Curtis basis
  std::unique_ptr<base::SBsplineModifiedClenshawCurtisBase> modBsplineClenshawCurtisBasis;
  /// fundamental spline basis
  std::unique_ptr<base::SFundamentalSplineBase> fundamentalSplineBasis;
  /// modified fundamental spline basis
  std::unique_ptr<base::SFundamentalSplineModifiedBase> modFundamentalSplineBasis;
  /// linear basis
  std::unique_ptr<base::SLinearBase> linearBasis;
  /// linear boundary basis
  std::unique_ptr<base::SLinearBoundaryBase> linearL0BoundaryBasis;
  /// linear Clenshaw-Curtis basis
  std::unique_ptr<base::SLinearClenshawCurtisBase> linearClenshawCurtisBasis;
  /// linear Clenshaw-Curtis boundary basis
  std::unique_ptr<base::SLinearClenshawCurtisBoundaryBase> linearClenshawCurtisBoundaryBasis;
  /// modified linear basis
  std::unique_ptr<base::SLinearModifiedBase> modLinearBasis;
  /// wavelet basis
  std::unique_ptr<base::SWaveletBase> waveletBasis;
  /// wavelet boundary basis
  std::unique_ptr<base::SWaveletBoundaryBase> waveletBoundaryBasis;
  /// modified wavelet basis
  std::unique_ptr<base::SWaveletModifiedBase> modWaveletBasis;
  /// not-a-knot B-spline Boundary basis
  std::unique_ptr<base::SNakBsplineBoundaryBase> nakBsplineBoundaryBasis;
  /// not-a-knot B-spline Boundary combigrid basis
  std::unique_ptr<base::SNakBsplineBoundaryCombigridBase> nakBsplineBoundaryCombigridBasis;
  /// not-a-knot B-spline Boundary combigrid basis
  std::unique_ptr<base::SNakBsplineModifiedBase> nakBsplineModifiedBasis;
  /// mod poly basis
  std::unique_ptr<base::SPolyModifiedBase> modPolyBasis;
  /// poly basis
  std::unique_ptr<base::SPolyBase> polyBasis;
  /// poly boundary basis
  std::unique_ptr<base::SPolyBoundaryBase> polyBoundaryBasis;
  /// nak Bspline basis
  std::unique_ptr<base::SNakBsplineBase> nakBsplineBasis;
  /// nak Bspline extended basis
  std::unique_ptr<base::SNakBsplineExtendedBase> nakBsplineExtendedBasis;

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
    NAK_BSPLINEBOUNDARY_COMBIGRID,
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
    } else if (basisType == NAK_BSPLINEBOUNDARY_COMBIGRID) {
      return evalNakBsplineBoundaryCombigridFunctionAtGridPoint(basisI, pointJ);
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
   * @return          value of the basisI-th not-a-knot B-spline combigrid basis function
   *                  at the pointJ-th grid point
   */
  inline double evalNakBsplineBoundaryCombigridFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = nakBsplineBoundaryCombigridBasis->eval(
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
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATIONSLE_HPP */
