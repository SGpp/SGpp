// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLagrangeNotAKnotSplineBoundary.hpp>

#include <sgpp/globaldef.hpp>

#include <deque>
#include <limits>
#include <vector>

namespace sgpp {
namespace base {

namespace {

void discoverGrid(size_t dim, size_t N, GridStorage::grid_iterator& index,
                  GridStorage& storage1D, std::vector<size_t>& k1D) {
  const size_t k = index.seq();

  if (k >= N) {
    return;
  }

  level_t l;
  index_t i;
  GridPoint point1D(1);

  index.get(dim, l, i);
  point1D.set(0, l, i);
  storage1D.insert(point1D);
  k1D.push_back(k);

  if (l >= 1) {
    index.set(dim, l + 1, 2 * i - 1);
    discoverGrid(dim, N, index, storage1D, k1D);

    index.set(dim, l + 1, 2 * i + 1);
    discoverGrid(dim, N, index, storage1D, k1D);
  }

  index.set(dim, l, i);
}

inline double hermiteInterpolation(double x, double a, double b,
                                   double fa, double da, double fb, double db) {
    x = (x - a) / (b - a);
    da *= b - a;
    db *= b - a;

    const double x2 = x * x;
    const double x3 = x2 * x;
    const double h0 = 2.0 * x3 - 3.0 * x2 + 1.0;
    const double h1 = x3 - 2.0 * x2 + x;
    const double h2 = x3 - x2;
    const double h3 = -2.0 * x3 + 3.0 * x2;

    return fa * h0 + da * h1 + db * h2 + fb * h3;
}

inline double hermiteInterpolationDeriv(double x, double a, double b,
                                        double fa, double da, double fb, double db) {
    x = (x - a) / (b - a);
    fa /= b - a;
    fb /= b - a;

    const double x2 = x * x;
    const double h0 = 6.0 * x2 - 6.0 * x;
    const double h1 = 3.0 * x2 - 4.0 * x + 1.0;
    const double h2 = 3.0 * x2 - 2.0 * x;
    const double h3 = -6.0 * x2 + 6.0 * x;

    return fa * h0 + da * h1 + db * h2 + fb * h3;
}

inline void setPointLevelEvenIndex(GridPoint& point,
                                   level_t l, index_t i) {
    level_t l2 = l;
    index_t i2 = i;

    while ((l2 > 0) && (i2 % 2 == 0)) {
        l2--;
        i2 /= 2;
    }

    point.set(0, l2, i2);
}

inline void hh1DLevelZero(const GridStorage& gridStorage,
                          const DataVector& fX,
                          DataVector& alpha,
                          DataVector& dX,
                          const size_t& N,
                          GridPoint& point,
                          std::deque<size_t>& kListNext) {
    point.getLeftLevelZero(0);
    const size_t kl = gridStorage.getSequenceNumber(point);
    point.getRightLevelZero(0);
    const size_t kr = gridStorage.getSequenceNumber(point);

    const double fl = fX[kl];
    const double fr = fX[kr];

    alpha[kl] = fl;
    alpha[kr] = fr;
    dX[kl] = fr - fl;
    dX[kr] = fr - fl;

    {
        point.getRoot(0);
        const size_t kChild = gridStorage.getSequenceNumber(point);

        if (kChild < N) {
            kListNext.push_back(kChild);
        }
    }
}

inline void hh1DForwardSubstitution(const GridStorage& gridStorage,
                                    const DataVector& fX,
                                    DataVector& alpha,
                                    DataVector& dX,
                                    const size_t& N,
                                    GridPoint& point,
                                    const std::deque<size_t>& kListCurrent,
                                    std::deque<size_t>& kListNext,
                                    const level_t& l,
                                    const index_t& hInv,
                                    std::vector<double>& matrixLeft,
                                    std::vector<double>& matrixCenter,
                                    std::vector<double>& matrixRight) {
  const double NaN = std::numeric_limits<double>::quiet_NaN();

  for (size_t q = 0; q < kListCurrent.size(); q++) {
    const size_t k = kListCurrent[q];
    point = gridStorage[k];
    const index_t i = point.getIndex(0);
    const double x = gridStorage.getCoordinate(point, 0);

    if (l == 1) {
      matrixLeft[q]   = NaN;
      matrixCenter[q] = 1.0;
      matrixRight[q]  = NaN;
    } else if ((l == 2) && (i == 1)) {
      matrixLeft[q]   = NaN;
      matrixCenter[q] = 9.0/26.0;
      matrixRight[q]  = -9.0/104.0;
    } else if ((l == 2) && (i == 3)) {
      matrixLeft[q]   = -9.0/104.0;
      matrixCenter[q] = 9.0/26.0;
      matrixRight[q]  = NaN;
    } else if (i == 1) {
      matrixLeft[q]   = NaN;
      matrixCenter[q] = 9.0/28.0;
      matrixRight[q]  = -1.0/56.0;
    } else if (i == 3) {
      matrixLeft[q]   = -3.0/16.0;
      matrixCenter[q] = 19.0/32.0;
      matrixRight[q]  = -1.0/24.0;
    } else if (i + 3 == hInv) {
      matrixLeft[q]   = -1.0/24.0;
      matrixCenter[q] = 19.0/32.0;
      matrixRight[q]  = -3.0/16.0;
    } else if (i + 1 == hInv) {
      matrixLeft[q]   = -1.0/56.0;
      matrixCenter[q] = 9.0/28.0;
      matrixRight[q]  = NaN;
    } else {
      matrixLeft[q]   = -1.0/24.0;
      matrixCenter[q] = 7.0/12.0;
      matrixRight[q]  = -1.0/24.0;
    }

    {
      setPointLevelEvenIndex(point, l, i - 1);
      const size_t kl = gridStorage.getSequenceNumber(point);
      const double xl = gridStorage.getCoordinate(point, 0);

      setPointLevelEvenIndex(point, l, i + 1);
      const size_t kr = gridStorage.getSequenceNumber(point);
      const double xr = gridStorage.getCoordinate(point, 0);

      const double ftx = hermiteInterpolation(x, xl, xr, fX[kl], dX[kl],
                                              fX[kr], dX[kr]);
      dX[k] = hermiteInterpolationDeriv(x, xl, xr, fX[kl], dX[kl], fX[kr], dX[kr]);
      alpha[k] = fX[k] - ftx;
    }

    if (q > 0) {
      const size_t ql = q - 1;
      const size_t kl = kListCurrent[ql];

      if (gridStorage[kl].getIndex(0) == i - 2) {
        const double rowFactor = matrixRight[ql] / matrixCenter[ql];
        alpha[k] -= rowFactor * alpha[kl];
        matrixCenter[q] -= rowFactor * matrixLeft[q];
        matrixRight[ql] = 0.0;
      }
    }

    {
      point.set(0, l + 1, 2 * i - 1);
      const size_t kChild = gridStorage.getSequenceNumber(point);

      if (kChild < N) {
        kListNext.push_back(kChild);
      }
    }

    {
      point.set(0, l + 1, 2 * i + 1);
      const size_t kChild = gridStorage.getSequenceNumber(point);

      if (kChild < N) {
        kListNext.push_back(kChild);
      }
    }

    point.set(0, l, i);
  }
}

inline void hh1DBackwardSubstitution(const GridStorage& gridStorage,
                                     DataVector& alpha,
                                     const std::deque<size_t>& kListCurrent,
                                     std::vector<double>& matrixLeft,
                                     std::vector<double>& matrixCenter) {
  for (size_t q = kListCurrent.size(); q-- > 0; ) {
    const size_t k = kListCurrent[q];
    const index_t i = gridStorage[k].getIndex(0);

    if (q < kListCurrent.size() - 1) {
      const size_t qr = q + 1;
      const size_t kr = kListCurrent[qr];

      if (gridStorage[kr].getIndex(0) == i + 2) {
        const double rowFactor = matrixLeft[qr] / matrixCenter[qr];
        alpha[k] -= rowFactor * alpha[kr];
        matrixLeft[qr] = 0.0;
      }
    }

    alpha[k] /= matrixCenter[q];
    matrixCenter[q] = 1.0;
  }
}

inline void hh1DUpdateDerivatives(const GridStorage& gridStorage,
                                  const DataVector& alpha,
                                  DataVector& dX,
                                  const size_t& N,
                                  GridPoint& point,
                                  const std::deque<size_t>& kListCurrent,
                                  const level_t& l,
                                  const index_t& hInv,
                                  std::vector<double>& derivValues) {
  const double NaN = std::numeric_limits<double>::quiet_NaN();

  for (size_t q = 0; q < kListCurrent.size(); q++) {
    // k = basis function, k2 = grid point
    const size_t k = kListCurrent[q];
    point = gridStorage[k];
    const index_t i = point.getIndex(0);

    const double innerDeriv = static_cast<double>(hInv);

    if (l == 1) {
      derivValues[0] = NaN;
      derivValues[1] = NaN;
      derivValues[2] = 4.0;
      derivValues[3] = 0.0;
      derivValues[4] = -4.0;
      derivValues[5] = NaN;
      derivValues[6] = NaN;
    } else if ((l == 2) && (i == 1)) {
      derivValues[0] = NaN;
      derivValues[1] = NaN;
      derivValues[2] = 57.0/13.0;
      derivValues[3] = -21.0/26.0;
      derivValues[4] = -15.0/13.0;
      derivValues[5] = 3.0/13.0;
      derivValues[6] = 3.0/13.0;
    } else if ((l == 2) && (i == 3)) {
      derivValues[0] = -3.0/13.0;
      derivValues[1] = -3.0/13.0;
      derivValues[2] = 15.0/13.0;
      derivValues[3] = 21.0/26.0;
      derivValues[4] = -57.0/13.0;
      derivValues[5] = NaN;
      derivValues[6] = NaN;
    } else if (i == 1) {
      derivValues[0] = NaN;
      derivValues[1] = NaN;
      derivValues[2] = 15.0/14.0 * innerDeriv;
      derivValues[3] = -3.0/14.0 * innerDeriv;
      derivValues[4] = -3.0/14.0 * innerDeriv;
      derivValues[5] = 3.0/56.0 * innerDeriv;
      derivValues[6] = 0.0;
    } else if (i == 3) {
      derivValues[0] = -1.0/8.0 * innerDeriv;
      derivValues[1] = -1.0/8.0 * innerDeriv;
      derivValues[2] = 5.0/8.0 * innerDeriv;
      derivValues[3] = -1.0/32.0 * innerDeriv;
      derivValues[4] = -1.0/2.0 * innerDeriv;
      derivValues[5] = 1.0/8.0 * innerDeriv;
      derivValues[6] = 0.0;
    } else if (i + 3 == hInv) {
      derivValues[0] = 0.0;
      derivValues[1] = -1.0/8.0 * innerDeriv;
      derivValues[2] = 1.0/2.0 * innerDeriv;
      derivValues[3] = 1.0/32.0 * innerDeriv;
      derivValues[4] = -5.0/8.0 * innerDeriv;
      derivValues[5] = 1.0/8.0 * innerDeriv;
      derivValues[6] = 1.0/8.0 * innerDeriv;
    } else if (i + 1 == hInv) {
      derivValues[0] = 0.0;
      derivValues[1] = -3.0/56.0 * innerDeriv;
      derivValues[2] = 3.0/14.0 * innerDeriv;
      derivValues[3] = 3.0/14.0 * innerDeriv;
      derivValues[4] = -15.0/14.0 * innerDeriv;
      derivValues[5] = NaN;
      derivValues[6] = NaN;
    } else {
      derivValues[0] = 0.0;
      derivValues[1] = -1.0/8.0 * innerDeriv;
      derivValues[2] = 1.0/2.0 * innerDeriv;
      derivValues[3] = 0.0;
      derivValues[4] = -1.0/2.0 * innerDeriv;
      derivValues[5] = 1.0/8.0 * innerDeriv;
      derivValues[6] = 0.0;
    }

    for (index_t j = 0; j <= 3; j++) {
      if (i >= j) {
        setPointLevelEvenIndex(point, l, i - j);
        const size_t k2 = gridStorage.getSequenceNumber(point);

        if (k2 < N) {
            dX[k2] += alpha[k] * derivValues[3 - j];
        }
      }

      if ((j > 0) && (i + j <= hInv)) {
        setPointLevelEvenIndex(point, l, i + j);
        const size_t k2 = gridStorage.getSequenceNumber(point);

        if (k2 < N) {
            dX[k2] += alpha[k] * derivValues[3 + j];
        }
      }
    }
  }
}

inline void linearHierarchize1D(GridStorage& gridStorage,
                                const DataVector& fX,
                                DataVector& alpha) {
  const size_t N = gridStorage.getSize();
  GridPoint point(1);
  alpha.resize(N);

  for (size_t k = 0; k < N; k++) {
    point = gridStorage[k];
    const level_t l = point.getLevel(0);
    const index_t i = point.getIndex(0);

    setPointLevelEvenIndex(point, l, i - 1);
    const size_t kl = gridStorage.getSequenceNumber(point);

    setPointLevelEvenIndex(point, l, i + 1);
    const size_t kr = gridStorage.getSequenceNumber(point);

    alpha[k] = fX[k] - (fX[kl] + fX[kr]) / 2.0;

    point.set(0, l, i);
  }
}

inline void hermiteHierarchize1D(GridStorage& gridStorage,
                                 const DataVector& fX,
                                 DataVector& alpha) {
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const size_t N = gridStorage.getSize();
  alpha.resize(N);
  alpha.setAll(NaN);
  DataVector dX(N, NaN);
  GridPoint point(1);
  std::deque<size_t> kListCurrent;
  std::deque<size_t> kListNext;
  const size_t lMax = gridStorage.getMaxLevel();
  std::vector<double> derivValues(7);

  for (level_t l = 0; l <= lMax; l++) {
    const index_t hInv = static_cast<index_t>(1) << l;

    kListCurrent = kListNext;
    kListNext.clear();

    std::vector<double> matrixLeft(kListCurrent.size());
    std::vector<double> matrixCenter(kListCurrent.size());
    std::vector<double> matrixRight(kListCurrent.size());

    if (l == 0) {
      hh1DLevelZero(gridStorage, fX, alpha, dX, N, point, kListNext);
    } else {
      hh1DForwardSubstitution(gridStorage, fX, alpha, dX, N, point, kListCurrent, kListNext,
                              l, hInv, matrixLeft, matrixCenter, matrixRight);
      hh1DBackwardSubstitution(gridStorage, alpha, kListCurrent, matrixLeft, matrixCenter);
      hh1DUpdateDerivatives(gridStorage, alpha, dX, N, point, kListCurrent,
                            l, hInv, derivValues);
    }
  }
}

}  // end namespace

HierarchisationLagrangeNotAKnotSplineBoundary::HierarchisationLagrangeNotAKnotSplineBoundary(
  GridStorage& storage) : storage(storage) {
}

HierarchisationLagrangeNotAKnotSplineBoundary::~HierarchisationLagrangeNotAKnotSplineBoundary() {
}

void HierarchisationLagrangeNotAKnotSplineBoundary::operator()(DataVector& source,
    DataVector& result, grid_iterator& index, size_t dim) {
  const size_t N = storage.getSize();
  GridStorage storage1D(1);
  std::vector<size_t> k1D;

  index.resetToLeftLevelZero(dim);
  discoverGrid(dim, N, index, storage1D, k1D);
  index.resetToRightLevelZero(dim);
  discoverGrid(dim, N, index, storage1D, k1D);
  index.resetToLevelOne(dim);
  discoverGrid(dim, N, index, storage1D, k1D);

  const size_t N1D = storage1D.getSize();
  DataVector source1D(N1D);
  DataVector result1D(N1D);

  for (size_t q = 0; q < k1D.size(); q++) {
    source1D[q] = source[k1D[q]];
    result1D[q] = result[k1D[q]];
  }

  hermiteHierarchize1D(storage1D, source1D, result1D);

  for (size_t q = 0; q < k1D.size(); q++) {
    result[k1D[q]] = result1D[q];
  }
}

}  // namespace base
}  // namespace sgpp
