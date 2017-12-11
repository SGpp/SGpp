// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>

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

// ToDo (rehmemk) Test other strategies for outer points for example placing them uniform using the
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
  // On levels 0,1 and 2 only Lagrange polynomials and not nak B splines are used. no need for extra
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
  sgpp::combigrid::GridFunction gf([=](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
    size_t numDimensions = grid->getDimension();
    // stores the values of the objective function
    auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);

    // ToDo (rehmemk) B-spline interpolation and B-spline quadrature can be mixed (one in one
    // dimension the other in another dimension and so on). Combining B-splines and other basis
    // functions hast not been tested yet.

    sgpp::combigrid::CombiEvaluators::Collection interpolEvaluators(
        numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));

    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double>>(numDimensions);
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

      std::vector<std::vector<double>> basisValues;
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
