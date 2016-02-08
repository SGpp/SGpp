// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <string>

namespace SGPP {
namespace datadriven {

DataMiningConfiguration::DataMiningConfiguration(): json::JSON() {
}

DataMiningConfiguration::DataMiningConfiguration(const std::string& fileName):
  json::JSON(fileName) {
}

DataMiningConfiguration* DataMiningConfiguration::clone() {
	DataMiningConfiguration* clone = new DataMiningConfiguration(*this);
  return clone;
}

SGPP::base::GridType DataMiningConfiguration::stringToGridType(std::string& gridType) {
  if (gridType.compare("Linear") == 0) {
    return SGPP::base::GridType::Linear;
  } else if (gridType.compare("LinearStretched") == 0) {
    return SGPP::base::GridType::LinearStretched;
  } else if (gridType.compare("LinearL0Boundary") == 0) {
    return SGPP::base::GridType::LinearL0Boundary;
  } else if (gridType.compare("LinearBoundary") == 0) {
    return SGPP::base::GridType::LinearBoundary;
  } else if (gridType.compare("LinearStretchedBoundary") == 0) {
    return SGPP::base::GridType::LinearStretchedBoundary;
  } else if (gridType.compare("LinearTruncatedBoundary") == 0) {
    return SGPP::base::GridType::LinearTruncatedBoundary;
  } else if (gridType.compare("ModLinear") == 0) {
    return SGPP::base::GridType::ModLinear;
  } else if (gridType.compare("Poly") == 0) {
    return SGPP::base::GridType::Poly;
  } else if (gridType.compare("PolyBoundary") == 0) {
    return SGPP::base::GridType::PolyBoundary;
  } else if (gridType.compare("ModPoly") == 0) {
    return SGPP::base::GridType::ModPoly;
  } else if (gridType.compare("ModWavelet") == 0) {
    return SGPP::base::GridType::ModWavelet;
  } else if (gridType.compare("ModBspline") == 0) {
    return SGPP::base::GridType::ModBspline;
  } else if (gridType.compare("Prewavelet") == 0) {
    return SGPP::base::GridType::Prewavelet;
  } else if (gridType.compare("SquareRoot") == 0) {
    return SGPP::base::GridType::SquareRoot;
  } else if (gridType.compare("Periodic") == 0) {
    return SGPP::base::GridType::Periodic;
  } else if (gridType.compare("LinearClenshawCurtis") == 0) {
    return SGPP::base::GridType::LinearClenshawCurtis;
  } else if (gridType.compare("Bspline") == 0) {
    return SGPP::base::GridType::Bspline;
  } else if (gridType.compare("BsplineBoundary") == 0) {
    return SGPP::base::GridType::BsplineBoundary;
  } else if (gridType.compare("BsplineClenshawCurtis") == 0) {
    return SGPP::base::GridType::BsplineClenshawCurtis;
  } else if (gridType.compare("Wavelet") == 0) {
    return SGPP::base::GridType::Wavelet;
  } else if (gridType.compare("WaveletBoundary") == 0) {
    return SGPP::base::GridType::WaveletBoundary;
  } else if (gridType.compare("FundamentalSpline") == 0) {
    return SGPP::base::GridType::FundamentalSpline;
  } else if (gridType.compare("ModFundamentalSpline") == 0) {
    return SGPP::base::GridType::ModFundamentalSpline;
  } else if (gridType.compare("ModBsplineClenshawCurtis") == 0) {
    return SGPP::base::GridType::ModBsplineClenshawCurtis;
  } else if (gridType.compare("LinearStencil") == 0) {
    return SGPP::base::GridType::LinearStencil;
  } else if (gridType.compare("ModLinearStencil") == 0) {
    return SGPP::base::GridType::ModLinearStencil;
  } else {
    throw SGPP::base::application_exception("grid type is unknown");
  }
}

SGPP::datadriven::RegularizationType DataMiningConfiguration::stringToRegularizationType(std::string& regularizationType) {
  if (regularizationType.compare("Identity") == 0) {
    return SGPP::datadriven::RegularizationType::Identity;
  } else if (regularizationType.compare("Laplace") == 0) {
    return SGPP::datadriven::RegularizationType::Laplace;
  } else {
    throw SGPP::base::application_exception("regularization type is unknown");
  }
}

SGPP::solver::SLESolverType DataMiningConfiguration::stringToSolverType(std::string& solverType) {
  if (solverType.compare("CG")) {
    return SGPP::solver::SLESolverType::CG;
  } else if (solverType.compare("BiCGSTAB")) {
    return SGPP::solver::SLESolverType::BiCGSTAB;
  } else {
    throw SGPP::base::application_exception("solver type is unknown");
  }
}

}  // namespace base
}  // namespace SGPP

