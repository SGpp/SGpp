// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>

// ToDo (rehmemk) Test other strategies for outer points for example placing them uniform using
// the
// distance of the last point inside to the boundary

std::vector<double> createdeg1Knots(std::vector<double> const& xValues) {
  size_t degree = 1;

  // this offset works only for odd degrees. if even degrees shall be supported it must be
  // generalized
  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * offset, 0);
  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset; i++) {
    xi[offset - i - 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }
  return xi;
}

std::vector<double> createdeg3NakKnots(std::vector<double> const& xValues) {
  size_t degree = 3;
  // On levels 0 and 1 only Lagrange polynomials and not nak B splines are used. no need for extra
  // points
  if (xValues.size() < 5) {
    return xValues;
  }
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * offset + 2, 0);

  xi.insert(xi.begin() + offset + 1, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 1; i++) {
    xi[offset - i] = -xValues.at(i + 1) + 2 * xValues[0];
    xi[xValues.size() + offset + i + 1] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues.at(xValues.size() - i - 2));
  }

  xi.erase(xi.begin() + offset + 2);
  xi.erase(xi.end() - offset - 3);
  return xi;
}

std::vector<double> createdeg5NakKnots(std::vector<double> const& xValues) {
  size_t degree = 5;
  // On levels 0,1 and 2 only Lagrange polynomials and not nak B splines are used. no need for
  // extra
  // points
  if (xValues.size() < 9) {
    return xValues;
  }
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * (offset + 2), 0);

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
  return xi;
}

// creates the not a knot knots xi from the 1D grid points xValues.
// xi then is used for the BsplineQuadrature and Bspline ScalarProduct to identify the supports of
// the nakBsplines
std::vector<double> createNakKnots(std::vector<double> const& xValues, size_t const& degree) {
  if (degree == 1) {
    return createdeg1Knots(xValues);
  } else if (degree == 3) {
    return createdeg3NakKnots(xValues);
  } else if (degree == 5) {
    return createdeg5NakKnots(xValues);
  } else {
    throw std::invalid_argument("BSplineRoutines: only B-Spline degrees 1,3 and 5 supported");
  }
}

// ToDo (rehmemk) use unidirectional principle instead of global SLE solving to speed this up?!

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

    sgpp::base::FullSLE sle(A);

    sgpp::base::sle_solver::Auto solver;
    sgpp::base::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

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

size_t getUniqueIndex(size_t level, size_t index) {
  size_t offset = 0;
  if (level == 0) {
    offset = 0;
  } else if (level == 1) {
    offset = 1;
  } else if (level == 2) {
    offset = 4;
  } else {
    offset = 4 + static_cast<size_t>(std::pow(2, level - 1));
  }
  return offset + index;
}

sgpp::combigrid::GridFunction BSplineTensorCoefficientGridFunction(
    sgpp::combigrid::MultiFunction func, sgpp::combigrid::CombiHierarchies::Collection grids,
    size_t degree) {
  sgpp::combigrid::GridFunction gf(
      [func, grids](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
    size_t numDimensions = grid->getDimension();
    // stores the values of the objective function
    auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);

    std::cout << "-----------------------------------------------" << std::endl;
    sgpp::combigrid::CombiEvaluators::TensorCollection interpolEvaluators(
        numDimensions, sgpp::combigrid::CombiEvaluators::tensorBSplineInterpolation());

    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double>>(numDimensions);
    auto level = grid->getLevel();
    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }

    std::vector<bool> orderingConfiguration;

    sgpp::combigrid::CombiEvaluators::TensorCollection evalCopy(numDimensions);
    for (size_t dim = 0; dim < numDimensions; ++dim) {
      evalCopy[dim] = interpolEvaluators[dim]->cloneLinear();
      bool needsSorted = evalCopy[dim]->needsOrderedPoints();
      auto gridPoints = grids[dim]->getPoints(level[dim], needsSorted);
      orderingConfiguration.push_back(needsSorted);
      evalCopy[dim]->setGridPoints(gridPoints);
    }
    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
    A.setAll(0.0);
    sgpp::base::DataVector coefficients_sle(numGridPoints);
    sgpp::base::DataVector functionValues(numGridPoints);

    // Creates an iterator that yields the multi-indices of all grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());

    auto funcIter = funcStorage->getGuidedIterator(level, it, orderingConfiguration);

    // assemble interpolation matrix from tensors and rhs
    for (size_t ixEvalPoints = 0; funcIter->isValid(); ++ixEvalPoints, funcIter->moveToNext()) {
      // assemble row ixEvalPoints of interpolation matrix
      // load function values
      functionValues[ixEvalPoints] = funcIter->value();

      // load basis evaluations per dimension from tensor
      sgpp::combigrid::MultiIndex ix = funcIter->getMultiIndex();
      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
      for (size_t ixBasisFunctions = 0; innerIter.isValid();
           ++ixBasisFunctions, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t idim = 0; idim < numDimensions; ++idim) {
          auto iy = sgpp::combigrid::MultiIndex{getUniqueIndex(level[idim], innerIndex[idim])};
          splineValue *= evalCopy[idim]->getBasisValues()[ix[idim]].get(iy).getValue();
        }
        A.set(ixBasisFunctions, ixEvalPoints, splineValue);
      }
    }

    sgpp::base::FullSLE sle(A);
    sgpp::base::sle_solver::Auto solver;
    sgpp::base::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

    if (!solved) {
      throw sgpp::base::application_exception(
          "BSplineRoutines::BSplineTensorCoefficientGridFunction - interpolation matrix is "
          "singular.");
    }

    //    std::cout << "coeff: " << std::endl;
    //    std::cout << "[";
    //    for (size_t i = 0; i < coefficients_sle.size() - 1; i++) {
    //      std::cout << std::setw(15) << std::setprecision(10) << coefficients_sle[i] << ", ";
    //    }
    //    std::cout << coefficients_sle[coefficients_sle.size() - 1] << "]" << std::endl;
    //    std::cout << "----------------------------------" << std::endl;

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
