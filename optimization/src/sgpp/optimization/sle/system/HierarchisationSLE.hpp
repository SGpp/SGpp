// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATIONSLE_HPP
#define SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATIONSLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>

#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <memory>
#include <functional>

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
  explicit HierarchisationSLE(base::Grid& grid) : HierarchisationSLE(grid, grid.getStorage()) {}

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
      : CloneableSLE(), grid(grid), gridStorage(gridStorage), basis(nullptr),
        p(0), pp1h(0), pp1hDbl(0.0) {
    if (grid.getType() == base::GridType::Bspline) {
      p = dynamic_cast<base::BsplineGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SBsplineBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x > h * (static_cast<double>(i) - pp1hDbl)) &&
               (x < h * (static_cast<double>(i) + pp1hDbl));
      };
    } else if (grid.getType() == base::GridType::BsplineBoundary) {
      p = dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SBsplineBoundaryBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x > h * (static_cast<double>(i) - pp1hDbl)) &&
               (x < h * (static_cast<double>(i) + pp1hDbl));
      };
    } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
      p = dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SBsplineClenshawCurtisBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const base::index_t hInv = static_cast<base::index_t>(1) << l;

        if (i >= pp1h) {
          const double xl = base::ClenshawCurtisTable::getInstance().getPoint(l, i - pp1h);

          if (x <= xl) {
            return false;
          }
        }

        if (i + pp1h <= hInv) {
          const double xr = base::ClenshawCurtisTable::getInstance().getPoint(l, i + pp1h);

          if (x >= xr) {
            return false;
          }
        }

        return true;
      };
    } else if (grid.getType() == base::GridType::ModBspline) {
      p = dynamic_cast<base::ModBsplineGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SBsplineModifiedBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x > h * (static_cast<double>(i) - pp1hDbl)) &&
               (x < h * (static_cast<double>(i) + pp1hDbl));
      };
    } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
      p = dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SBsplineModifiedClenshawCurtisBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        return true;
        const base::index_t hInv = static_cast<base::index_t>(1) << l;

        if (i >= pp1h) {
          const double xl = base::ClenshawCurtisTable::getInstance().getPoint(l, i - pp1h);

          if (x <= xl) {
            return false;
          }
        }

        if (i + pp1h <= hInv) {
          const double xr = base::ClenshawCurtisTable::getInstance().getPoint(l, i + pp1h);

          if (x >= xr) {
            return false;
          }
        }

        return true;
      };
    } else if (grid.getType() == base::GridType::FundamentalSpline) {
      p = dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SFundamentalSplineBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double hInv = static_cast<double>(1 << l);
        const double xoh = x * hInv;

        if ((xoh == static_cast<double>(static_cast<base::index_t>(xoh))) &&
            (xoh != static_cast<double>(i))) {
          return false;
        } else {
          return true;
        }
      };
    } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
      p = dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree();
      basis = std::unique_ptr<base::SBasis>(new base::SFundamentalSplineModifiedBase(p));
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const base::index_t hInv = static_cast<base::index_t>(1) << l;
        const double hInvDbl = static_cast<double>(hInv);
        const double xoh = x * hInvDbl;

        if ((xoh == static_cast<double>(static_cast<base::index_t>(xoh))) &&
            (xoh != static_cast<double>(i)) &&
            !((i == 1) && (x == 0)) && !((i == hInv - 1) && (x == 1))) {
          return false;
        } else {
          return true;
        }
      };
    } else if (grid.getType() == base::GridType::Linear) {
      basis = std::unique_ptr<base::SBasis>(new base::SLinearBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x > h * (static_cast<double>(i) - 1.0)) &&
               (x < h * (static_cast<double>(i) + 1.0));
      };
    } else if (grid.getType() == base::GridType::LinearBoundary) {
      basis = std::unique_ptr<base::SBasis>(new base::SLinearBoundaryBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x > h * (static_cast<double>(i) - 1.0)) &&
               (x < h * (static_cast<double>(i) + 1.0));
      };
    } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
      basis = std::unique_ptr<base::SBasis>(new base::SLinearClenshawCurtisBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const base::index_t hInv = static_cast<base::index_t>(1) << l;

        if (i >= 1) {
          const double xl = base::ClenshawCurtisTable::getInstance().getPoint(l, i - 1);

          if (x <= xl) {
            return false;
          }
        }

        if (i + 1 <= hInv) {
          const double xr = base::ClenshawCurtisTable::getInstance().getPoint(l, i + 1);

          if (x >= xr) {
            return false;
          }
        }

        return true;
      };
    } else if (grid.getType() == base::GridType::ModLinear) {
      basis = std::unique_ptr<base::SBasis>(new base::SLinearModifiedBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x > h * (static_cast<double>(i) - 1.0)) &&
               (x < h * (static_cast<double>(i) + 1.0));
      };
    } else if (grid.getType() == base::GridType::Wavelet) {
      basis = std::unique_ptr<base::SBasis>(new base::SWaveletBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x >= h * (static_cast<double>(i) - 2.0)) &&
               (x <= h * (static_cast<double>(i) + 2.0));
      };
    } else if (grid.getType() == base::GridType::WaveletBoundary) {
      basis = std::unique_ptr<base::SBasis>(new base::SWaveletBoundaryBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x >= h * (static_cast<double>(i) - 2.0)) &&
               (x <= h * (static_cast<double>(i) + 2.0));
      };
    } else if (grid.getType() == base::GridType::ModWavelet) {
      basis = std::unique_ptr<base::SBasis>(new base::SWaveletModifiedBase());
      isPointIn1DSupportFunction = [this](base::level_t l, base::index_t i, double x) {
        const double h = 1.0 / static_cast<double>(1 << l);
        return (x >= h * (static_cast<double>(i) - 2.0)) &&
               (x <= h * (static_cast<double>(i) + 2.0));
      };
    } else {
      throw std::invalid_argument("Grid type not supported.");
    }

    pp1h = (static_cast<base::index_t>(p) + 1) / 2;
    pp1hDbl = static_cast<double>(pp1h);
  }

  /**
   * @param i     row index
   * @param j     column index
   * @return      whether the i-th grid point lies in the support of
   *              the j-th basis function
   */
  inline bool isMatrixEntryNonZero(size_t i, size_t j) override {
    return isGridPointInBasisFunctionSupport(j, i);
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
  typedef std::function<bool(base::level_t, base::index_t, double)> IsPointIn1DSupportFunction;

  /// sparse grid
  base::Grid& grid;
  /// grid storage
  base::GridStorage& gridStorage;
  /// basis
  std::unique_ptr<base::SBasis> basis;
  /// B-spline degree
  size_t p;
  /// (p + 1) / 2
  base::index_t pp1h;
  /// (p + 1) / 2
  double pp1hDbl;
  /// function to check if 1D point is in 1D basis support
  IsPointIn1DSupportFunction isPointIn1DSupportFunction;

  /**
   * @param basisI    basis function index
   * @param pointJ    grid point index
   * @return          value of the basisI-th basis function at the
   *                  pointJ-th grid point
   */
  inline double evalBasisFunctionAtGridPoint(size_t basisI, size_t pointJ) {
    const base::GridPoint& gpBasis = gridStorage[basisI];
    const base::GridPoint& gpPoint = gridStorage[pointJ];
    double result = 1.0;

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      const double result1d = basis->eval(
          gpBasis.getLevel(t), gpBasis.getIndex(t), gridStorage.getCoordinate(gpPoint, t));

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
   * @return          if false, the basis function is guaranteed to be zero at the grid point
   *                  (should be non-zero if true, but that's not guaranteed)
   */
  inline bool isGridPointInBasisFunctionSupport(size_t basisI, size_t pointJ) {
    const base::GridPoint& gpBasis = gridStorage[basisI];
    const base::GridPoint& gpPoint = gridStorage[pointJ];

    for (size_t t = 0; t < gridStorage.getDimension(); t++) {
      if (!isPointIn1DSupportFunction(gpBasis.getLevel(t), gpBasis.getIndex(t),
                                      gridStorage.getCoordinate(gpPoint, t))) {
        return false;
      }
    }

    return true;
  }
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATIONSLE_HPP */
