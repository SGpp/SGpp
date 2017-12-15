// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>

#include <vector>
#include <algorithm>

constexpr size_t log2(size_t n) { return ((n < 2) ? 1 : 1 + log2(n / 2)); }

// ToDo (rehmemk) has not been tested so far. Degree 5 still needs the level 0,1,2 changes
double expUniformNaKBspline(double const& x, size_t const& degree, size_t i,
                            std::vector<double> const& points) {
  size_t l = 0;
  if (points.size() > 1) {
    l = log2(points.size() - 2);
  }
  const size_t hInv = 1 << l;
  double t = x * static_cast<double>(hInv) - static_cast<double>(i);

  switch (degree) {
    case 1:
      return std::max(1.0 - std::abs(t), 0.0);

    // degree 3: polynomials on Level 0 and 1, nak Bsplines from Level 3 on
    case 3: {
      sgpp::base::BsplineBasis<size_t, size_t> bsplineBasis(3);
      if (l == 0) {
        if (i == 0) {
          // l = 0, i = 0
          std::cout << "one" << std::endl;
          return 1;
        }
        // 1 - x * x

      } else if (l == 1) {
        // Lagrange polynomials
        if (i == 0) {
          // l = 1, i = 0
          return 2 * x * x - 3 * x + 1;
        } else if (i == 1) {
          // l = 1, i = 1
          return 4 * x - 4 * x * x;
        } else {
          // l = 1, i = 2
          return 2 * x * x - x;
        }
      } else if ((i > 3) && (i < hInv - 3)) {
        // l >= 4, 3 < i < 2^l - 3
        return bsplineBasis.eval(l, i, x);
      } else {
        if (i > hInv / 2) {
          i = hInv - i;
          t *= -1.0;
        }

        if (i == 0) {
          // l >= 2, i = 0
          if ((t < 0.0) || (t > 2.0)) {
            return 0.0;
          } else {
            double result = -4.1666666666666664e-02;
            result = 2.5000000000000000e-01 + result * t;
            result = -5.0000000000000000e-01 + result * t;
            result = 3.3333333333333331e-01 + result * t;
            return result;
          }
        } else if ((l == 2) && (i == 1)) {
          // l = 2, i = 1
          if ((t < -1.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 1.0;
            double result = 1.0 / 10.0;
            result = -9.0 / 20.0 + result * t;
            result = 3.0 / 10.0 + result * t;
            result = 3.0 / 5.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 40.0;
            result = 3.0 / 20.0 + result * t;
            result = -3.0 / 10.0 + result * t;
            result = 1.0 / 5.0 + result * t;
            return result;
          }

        } else if (l == 2) {
          // l = 2, i = 2
          if ((t < -2.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 2.0;
            double result = -8.3333333333333329e-02;
            result = 2.0000000000000001e-01 + result * t;
            result = 2.0000000000000001e-01 + result * t;
            result = 6.6666666666666666e-02 + result * t;
            return result;
          } else {
            double result = 8.3333333333333329e-02;
            result = -2.9999999999999999e-01 + result * t;
            result *= t;
            result = 5.9999999999999998e-01 + result * t;
            return result;
          }
        } else if (i == 1) {
          // l >= 3, i = 1
          if ((t < -1.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 1.0;
            double result = 1.0 / 8.0;
            result = -1.0 / 2.0 + result * t;
            result = 1.0 / 4.0 + result * t;
            result = 7.0 / 12.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 12.0;
            result = 1.0 / 4.0 + result * t;
            result = -1.0 / 4.0 + result * t;
            result = 1.0 / 12.0 + result * t;
            return result;
          }
        } else if (i == 2) {
          // l >= 3, i = 2
          if ((t < -2.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 2.0;
            double result = -1.2500000000000000e-01;
            result = 2.5000000000000000e-01 + result * t;
            result = 2.5000000000000000e-01 + result * t;
            result = 8.3333333333333329e-02 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = 2.9166666666666669e-01;
            result = -5.0000000000000000e-01 + result * t;
            result = -2.5000000000000000e-01 + result * t;
            result = 5.8333333333333337e-01 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.2500000000000000e-01;
            result = 3.7500000000000000e-01 + result * t;
            result = -3.7500000000000000e-01 + result * t;
            result = 1.2500000000000000e-01 + result * t;
            return result;
          }
        } else {
          // l >= 3, i = 3
          if ((t < -3.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < -1.0) {
            t += 3.0;
            double result = 1.0 / 24.0;
            result *= t;
            result *= t;
            result *= t;
            return result;
          } else if (t < 0.0) {
            t += 1.0;
            double result = -3.0 / 8.0;
            result = 1.0 / 4.0 + result * t;
            result = 1.0 / 2.0 + result * t;
            result = 1.0 / 3.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = 11.0 / 24.0;
            result = -7.0 / 8.0 + result * t;
            result = -1.0 / 8.0 + result * t;
            result = 17.0 / 24.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 6.0;
            result = 1.0 / 2.0 + result * t;
            result = -1.0 / 2.0 + result * t;
            result = 1.0 / 6.0 + result * t;
            return result;
          }
        }
      }
    }
    // degree 5: Levels 0,1 and 2 polynomials, nak Bsplines from Level 3 on
    case 5: {
      sgpp::base::BsplineBasis<size_t, size_t> bsplineBasis(3);
      if (l == 0) {
        if (i == 0) {
          // l = 0, i = 0
          return 1 - x * x;  // 1 - x * x;
        } else {
          // l = 0, i = 1
          return x;
        }
      } else if (l == 1) {
        if (i == 1) {
          // l = 1, i = 1
          return 1;
        }
      } else if (l == 2) {
        if (i == 1) {
          // l = 2, i = 1 : cubic polynomial, 0 in 0,0.5,0.75 and 1 in 0.25
          return 32 * x * (x - 0.5) * (x - 0.75);
        } else if (i == 3) {
          // l = 2, i = 3 : quartic polynomial, 0 in 0,0.25,0.5,1 and 1 in 0.75
          return x * x * x * x;  // x * (x - 0.25) * (x - 0.5) * (x - 1) * (-128.0 / 3.0);
        }
      } else if ((i > 5) && (i < hInv - 5)) {
        // l >= 4, 5 < i < 2^l - 5
        return bsplineBasis.eval(l, i, x);
      } else {
        if (i > hInv / 2) {
          i = hInv - i;
          t *= -1.0;
        }

        if ((l == 3) && (i == 3)) {
          // l = 3, i = 3
          if ((t < -3.0) || (t > 5.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 3.0;
            double result = 107.0 / 30240.0;
            result = -17.0 / 756.0 + result * t;
            result = 1.0 / 378.0 + result * t;
            result = 37.0 / 378.0 + result * t;
            result = 109.0 / 756.0 + result * t;
            result = 253.0 / 3780.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -397.0 / 30240.0;
            result = 185.0 / 6048.0 + result * t;
            result = 155.0 / 3024.0 + result * t;
            result = -415.0 / 3024.0 + result * t;
            result = -1165.0 / 6048.0 + result * t;
            result = 2965.0 / 6048.0 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 233.0 / 30240.0;
            result = -53.0 / 1512.0 + result * t;
            result = 8.0 / 189.0 + result * t;
            result = 13.0 / 189.0 + result * t;
            result = -97.0 / 378.0 + result * t;
            result = 433.0 / 1890.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 4320.0;
            result = 1.0 / 288.0 + result * t;
            result = -1.0 / 48.0 + result * t;
            result = 1.0 / 16.0 + result * t;
            result = -3.0 / 32.0 + result * t;
            result = 9.0 / 160.0 + result * t;
            return result;
          }
        } else if (i == 1) {
          // l >= 3, i = 1
          if ((t < -1.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 2.0) {
            t += 1.0;
            double result = 1.0 / 504.0;
            result = -1.0 / 42.0 + result * t;
            result = 2.0 / 21.0 + result * t;
            result = -2.0 / 21.0 + result * t;
            result = -5.0 / 21.0 + result * t;
            result = 47.0 / 105.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 840.0;
            result = 1.0 / 168.0 + result * t;
            result = -1.0 / 84.0 + result * t;
            result = 1.0 / 84.0 + result * t;
            result = -1.0 / 168.0 + result * t;
            result = 1.0 / 840.0 + result * t;
            return result;
          }
        } else if (i == 3) {
          // l >= 4, i = 3
          if ((t < -3.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 3.0;
            double result = 1.0 / 252.0;
            result = -1.0 / 42.0 + result * t;
            result *= t;
            result = 2.0 / 21.0 + result * t;
            result = 1.0 / 7.0 + result * t;
            result = 1.0 / 15.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -23.0 / 1260.0;
            result = 1.0 / 28.0 + result * t;
            result = 1.0 / 14.0 + result * t;
            result = -5.0 / 42.0 + result * t;
            result = -1.0 / 4.0 + result * t;
            result = 163.0 / 420.0 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 19.0 / 1260.0;
            result = -1.0 / 18.0 + result * t;
            result = 2.0 / 63.0 + result * t;
            result = 8.0 / 63.0 + result * t;
            result = -2.0 / 9.0 + result * t;
            result = 34.0 / 315.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 252.0;
            result = 5.0 / 252.0 + result * t;
            result = -5.0 / 126.0 + result * t;
            result = 5.0 / 126.0 + result * t;
            result = -5.0 / 252.0 + result * t;
            result = 1.0 / 252.0 + result * t;
            return result;
          }
        } else {
          // l >= 4, i = 5
          if ((t < -5.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < -2.0) {
            t += 5.0;
            double result = 1.0 / 2520.0;
            result *= t;
            result *= t;
            result *= t;
            result *= t;
            result *= t;
            return result;
          } else if (t < -1.0) {
            t += 2.0;
            double result = -11.0 / 504.0;
            result = 1.0 / 168.0 + result * t;
            result = 1.0 / 28.0 + result * t;
            result = 3.0 / 28.0 + result * t;
            result = 9.0 / 56.0 + result * t;
            result = 27.0 / 280.0 + result * t;
            return result;
          } else if (t < 0.0) {
            t += 1.0;
            double result = 31.0 / 504.0;
            result = -13.0 / 126.0 + result * t;
            result = -10.0 / 63.0 + result * t;
            result = 2.0 / 63.0 + result * t;
            result = 25.0 / 63.0 + result * t;
            result = 121.0 / 315.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -181.0 / 2520.0;
            result = 103.0 / 504.0 + result * t;
            result = 11.0 / 252.0 + result * t;
            result = -113.0 / 252.0 + result * t;
            result = -61.0 / 504.0 + result * t;
            result = 1543.0 / 2520.0 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 11.0 / 280.0;
            result = -13.0 / 84.0 + result * t;
            result = 1.0 / 7.0 + result * t;
            result = 4.0 / 21.0 + result * t;
            result = -3.0 / 7.0 + result * t;
            result = 23.0 / 105.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 120.0;
            result = 1.0 / 24.0 + result * t;
            result = -1.0 / 12.0 + result * t;
            result = 1.0 / 12.0 + result * t;
            result = -1.0 / 24.0 + result * t;
            result = 1.0 / 120.0 + result * t;
            return result;
          }
        }
      }
    }

    default:
      return 0.0;
  }
}

double nonUniformBSpline(double const& x, size_t const& deg, size_t const& index,
                         std::vector<double> const& xi) {
  if (deg == 0) {
    // characteristic function of [xi[k], xi[k+1])
    return (((x >= xi[index]) && (x < xi[index + 1])) ? 1.0 : 0.0);
  } else if ((x < xi[index]) || (x >= xi[index + deg + 1])) {
    // out of support
    return 0.0;
  } else {
    // Cox-de-Boor recursion
    return (x - xi[index]) / (xi[index + deg] - xi[index]) *
               nonUniformBSpline(x, deg - 1, index, xi) +
           (1.0 - (x - xi[index + 1]) / (xi[index + deg + 1] - xi[index + 1])) *
               nonUniformBSpline(x, deg - 1, index + 1, xi);
  }
}

double LagrangePolynomial(double const& x, std::vector<double> const& xValues, size_t const& k) {
  double res = 1.0;
  for (size_t m = 0; m < xValues.size(); m++) {
    if (k != m) {
      res *= (x - xValues[m]) / (xValues[k] - xValues[m]);
    }
  }
  return res;
}

// ToDo (rehmemk) Test other strategies for outer points for example placing them uniform using
// the
// distance of the last point inside to the boundary

void createdeg1Knots(std::vector<double> const& xValues, std::vector<double>& xi) {
  size_t degree = 1;

  // this offset works only for odd degrees. if even degrees shall be supported it must be
  // generalized
  size_t offset = (degree + 1) / 2;
  xi.resize(2 * offset, 0);
  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset; i++) {
    xi[offset - i - 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }
}

void createdeg3NakKnots(std::vector<double> const& xValues, std::vector<double>& xi) {
  size_t degree = 3;
  // On levels 0 and 1 only Lagrange polynomials and not nak B splines are used. no need for extra
  // points
  if (xValues.size() < 5) {
    xi = xValues;
    return;
  }
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  xi.resize(2 * offset + 2, 0);

  xi.insert(xi.begin() + offset + 1, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 1; i++) {
    xi[offset - i] = -xValues.at(i + 1) + 2 * xValues[0];
    xi[xValues.size() + offset + i + 1] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues.at(xValues.size() - i - 2));
  }

  xi.erase(xi.begin() + offset + 2);
  xi.erase(xi.end() - offset - 3);
}

void createdeg5NakKnots(std::vector<double> const& xValues, std::vector<double>& xi) {
  size_t degree = 5;
  // On levels 0,1 and 2 only Lagrange polynomials and not nak B splines are used. no need for
  // extra
  // points
  if (xValues.size() < 9) {
    xi = xValues;
    return;
  }
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  xi.resize(2 * (offset + 2), 0);

  xi.insert(xi.begin() + offset + 2, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 2; i++) {
    xi[offset - i + 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i + 2] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }

  xi.erase(xi.begin() + offset + 3);
  xi.erase(xi.begin() + offset + 3);
  xi.erase(xi.end() - offset - 4);
  xi.erase(xi.end() - offset - 4);
}

void createNakKnots(std::vector<double> const& xValues, size_t const& degree,
                    std::vector<double>& xi) {
  if (degree == 1) {
    createdeg1Knots(xValues, xi);
  } else if (degree == 3) {
    createdeg3NakKnots(xValues, xi);
  } else if (degree == 5) {
    createdeg5NakKnots(xValues, xi);
  } else {
    throw std::invalid_argument("BSplineRoutines: only B-Spline degrees 1,3 and 5 supported");
  }
}

// ToDo (rehmemk) use unidirectional principle instead of global SLE solving to speed this up!

sgpp::combigrid::GridFunction BSplineCoefficientGridFunction(
    sgpp::combigrid::MultiFunction func, sgpp::combigrid::CombiHierarchies::Collection grids,
    size_t degree) {
  sgpp::combigrid::GridFunction gf([func, grids,
                                    degree](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
    size_t numDimensions = grid->getDimension();
    // stores the values of the objective function
    auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);

    // ToDo (rehmemk) B-spline interpolation and B-spline quadrature can be mixed (one in one
    // dimension the other in another dimension and so on). Combining B-splines and other basis
    // functions has not been tested yet.

    sgpp::combigrid::CombiEvaluators::Collection interpolEvaluators(
        numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));

    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double> >(numDimensions);
    auto level = grid->getLevel();
    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }

    std::vector<bool> orderingConfiguration;

    sgpp::combigrid::CombiEvaluators::Collection evalCopy(numDimensions);
    for (size_t dim = 0; dim < numDimensions; ++dim) {
      evalCopy[dim] = interpolEvaluators[dim]->cloneLinear();
      bool needsSorted = evalCopy[dim]->needsOrderedPoints();
      auto gridPoints = grids[dim]->getPoints(level[dim], needsSorted);
      orderingConfiguration.push_back(needsSorted);
      evalCopy[dim]->setGridPoints(gridPoints);
    }
    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
    sgpp::base::DataVector coefficients_sle(numGridPoints);
    sgpp::base::DataVector functionValues(numGridPoints);

    // Creates an iterator that yields the multi-indices of all grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());

    auto funcIter = funcStorage->getGuidedIterator(level, it, orderingConfiguration);

    for (size_t ixEvalPoints = 0; funcIter->isValid(); ++ixEvalPoints, funcIter->moveToNext()) {
      auto gridPoint = grid->getGridPoint(funcIter->getMultiIndex());
      functionValues[ixEvalPoints] = funcIter->value();

      std::vector<std::vector<double> > basisValues;
      for (size_t dim = 0; dim < numDimensions; ++dim) {
        evalCopy[dim]->setParameter(sgpp::combigrid::FloatScalarVector(gridPoint[dim]));
        auto basisValues1D = evalCopy[dim]->getBasisValues();
        // basis values at gridPoint
        std::vector<double> basisValues1D_vec(basisValues1D.size());
        for (size_t i = 0; i < basisValues1D.size(); i++) {
          basisValues1D_vec[i] = basisValues1D[i].value();
        }
        basisValues.push_back(basisValues1D_vec);
      }

      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
      for (size_t ixBasisFunctions = 0; innerIter.isValid();
           ++ixBasisFunctions, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t dim = 0; dim < numDimensions; ++dim) {
          splineValue *= basisValues[dim][innerIndex[dim]];
        }
        A.set(ixEvalPoints, ixBasisFunctions, splineValue);
      }
    }

    sgpp::optimization::FullSLE sle(A);
    sgpp::optimization::sle_solver::Auto solver;
    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

    /*
    std::cout << A.toString() << std::endl;
    std::cout << "fct: ";
    for (size_t i = 0; i < functionValues.size(); i++) {
      std::cout << functionValues[i] << " ";
    }
    std::cout << "\ncoeff: ";
    for (size_t i = 0; i < coefficients_sle.size(); i++) {
      std::cout << coefficients_sle[i] << " ";
    }
    std::cout << "\n";
    std::cout << "--------" << std::endl;
    */

    if (!solved) {
      exit(-1);
    }

    it.reset();
    for (size_t vecIndex = 0; it.isValid(); ++vecIndex, it.moveToNext()) {
      coefficientTree->set(it.getMultiIndex(), coefficients_sle[vecIndex]);
    }

    return coefficientTree;
  });
  return gf;
}
