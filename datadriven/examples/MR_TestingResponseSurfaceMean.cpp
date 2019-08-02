// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionNormal.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/datadriven/activeSubspaces/NakBsplineScalarProducts.hpp>
#include <sgpp/optimization/function/scalar/SparseGridResponseSurfaceBspline.hpp>

#include <iostream>
#include <random>

// double f(sgpp::base::DataVector v) { return sin(2 * M_PI * v[0]); }
double f(sgpp::base::DataVector v) { return sin(v[0]) * cos(v[1]); }

int main() {
  //  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  //  size_t degree = 3;
  //  size_t numDim = 2;
  //  size_t level = 8;
  //
  //  auto objectiveFunc = std::make_shared<sgpp::optimization::WrapperScalarFunction>(numDim, f);
  //  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineExtended;
  //  sgpp::base::DataVector lb = objectiveFunc->getLowerBounds();
  //  sgpp::base::DataVector ub = objectiveFunc->getUpperBounds();
  //  sgpp::optimization::SparseGridResponseSurfaceBspline reSurf(objectiveFunc, gridType, degree);
  //  reSurf.regular(level);
  //  //  auto pdf = std::make_shared<sgpp::datadriven::Uniform>(0, 1);
  //  auto pdf = std::make_shared<sgpp::base::DistributionNormal>(0.5, 0.1);
  //  size_t quadOrder = 20;
  //  double realMean = 0.416548687489679;
  //  double reSurfMean = reSurf.getMean(pdf, quadOrder);
  //
  //  //  sgpp::datadriven::Normal pdf_normal(0.5, 0.1);
  //  //  std::cout << pdf_normal.sample() << "\n";
  //
  //  std::cout << "#######################\n";
  //  sgpp::base::SNakBsplineExtendedBase base(degree);
  //  double basisMean = base.getMean(1, 1, pdf, quadOrder);
  //  double realBasisMean = 1.0;
  //
  //  std::cout << "\n\n";
  //  std::cout << "interpol error " << reSurf.l2Error(objectiveFunc, 10000) << "\n";
  //  std::cout << "mean " << reSurfMean << "\n";
  //  std::cout << "mean error " << reSurfMean - realMean << "\n";
  //
  //  //  std::cout << "basis mean " << basisMean << "\n";
  //  //  std::cout << "basis mean error " << realBasisMean - basisMean << "\n";

  return 0;
}
