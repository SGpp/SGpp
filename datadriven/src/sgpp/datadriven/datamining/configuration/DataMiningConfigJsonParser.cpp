// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigJsonParser.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

DataMiningConfigJsonParser::DataMiningConfigJsonParser() : json::JSON() {}

DataMiningConfigJsonParser::DataMiningConfigJsonParser(const std::string& fileName)
    : json::JSON(fileName) {}

DataMiningConfigJsonParser* DataMiningConfigJsonParser::clone() {
  DataMiningConfigJsonParser* clone = new DataMiningConfigJsonParser(*this);
  return clone;
}

sgpp::base::GridType DataMiningConfigJsonParser::stringToGridType(std::string& gridType) {
  if (gridType.compare("Linear") == 0) {
    return sgpp::base::GridType::Linear;
  } else if (gridType.compare("LinearStretched") == 0) {
    return sgpp::base::GridType::LinearStretched;
  } else if (gridType.compare("LinearL0Boundary") == 0) {
    return sgpp::base::GridType::LinearL0Boundary;
  } else if (gridType.compare("LinearBoundary") == 0) {
    return sgpp::base::GridType::LinearBoundary;
  } else if (gridType.compare("LinearStretchedBoundary") == 0) {
    return sgpp::base::GridType::LinearStretchedBoundary;
  } else if (gridType.compare("LinearTruncatedBoundary") == 0) {
    return sgpp::base::GridType::LinearTruncatedBoundary;
  } else if (gridType.compare("ModLinear") == 0) {
    return sgpp::base::GridType::ModLinear;
  } else if (gridType.compare("Poly") == 0) {
    return sgpp::base::GridType::Poly;
  } else if (gridType.compare("PolyBoundary") == 0) {
    return sgpp::base::GridType::PolyBoundary;
  } else if (gridType.compare("ModPoly") == 0) {
    return sgpp::base::GridType::ModPoly;
  } else if (gridType.compare("ModWavelet") == 0) {
    return sgpp::base::GridType::ModWavelet;
  } else if (gridType.compare("ModBspline") == 0) {
    return sgpp::base::GridType::ModBspline;
  } else if (gridType.compare("Prewavelet") == 0) {
    return sgpp::base::GridType::Prewavelet;
  } else if (gridType.compare("SquareRoot") == 0) {
    return sgpp::base::GridType::SquareRoot;
  } else if (gridType.compare("Periodic") == 0) {
    return sgpp::base::GridType::Periodic;
  } else if (gridType.compare("LinearClenshawCurtis") == 0) {
    return sgpp::base::GridType::LinearClenshawCurtis;
  } else if (gridType.compare("Bspline") == 0) {
    return sgpp::base::GridType::Bspline;
  } else if (gridType.compare("BsplineBoundary") == 0) {
    return sgpp::base::GridType::BsplineBoundary;
  } else if (gridType.compare("BsplineClenshawCurtis") == 0) {
    return sgpp::base::GridType::BsplineClenshawCurtis;
  } else if (gridType.compare("Wavelet") == 0) {
    return sgpp::base::GridType::Wavelet;
  } else if (gridType.compare("WaveletBoundary") == 0) {
    return sgpp::base::GridType::WaveletBoundary;
  } else if (gridType.compare("FundamentalSpline") == 0) {
    return sgpp::base::GridType::FundamentalSpline;
  } else if (gridType.compare("ModFundamentalSpline") == 0) {
    return sgpp::base::GridType::ModFundamentalSpline;
  } else if (gridType.compare("ModBsplineClenshawCurtis") == 0) {
    return sgpp::base::GridType::ModBsplineClenshawCurtis;
  } else if (gridType.compare("LinearStencil") == 0) {
    return sgpp::base::GridType::LinearStencil;
  } else if (gridType.compare("ModLinearStencil") == 0) {
    return sgpp::base::GridType::ModLinearStencil;
  } else {
    throw sgpp::base::application_exception("grid type is unknown");
  }
}

sgpp::datadriven::RegularizationType DataMiningConfigJsonParser::stringToRegularizationType(
    std::string& regularizationType) {
  if (regularizationType.compare("Identity") == 0) {
    return sgpp::datadriven::RegularizationType::Identity;
  } else if (regularizationType.compare("Laplace") == 0) {
    return sgpp::datadriven::RegularizationType::Laplace;
  } else {
    throw sgpp::base::application_exception("regularization type is unknown");
  }
}

sgpp::solver::SLESolverType DataMiningConfigJsonParser::stringToSolverType(
    std::string& solverType) {
  if (solverType.compare("CG")) {
    return sgpp::solver::SLESolverType::CG;
  } else if (solverType.compare("BiCGSTAB")) {
    return sgpp::solver::SLESolverType::BiCGSTAB;
  } else {
    throw sgpp::base::application_exception("solver type is unknown");
  }
}

}  // namespace datadriven
}  // namespace sgpp
