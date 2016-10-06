/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GridTypeParser.hpp
 *
 * Created on: Aug 18, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/base/exception/data_exception.hpp>

namespace sgpp {
namespace base {
class GridTypeParser {
 public:
  static GridType parse(const std::string& input) {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower.compare("linear") == 0) {
      return sgpp::base::GridType::Linear;
    } else if (inputLower.compare("linearstretched") == 0) {
      return sgpp::base::GridType::LinearStretched;
    } else if (inputLower.compare("linearl0boundary") == 0) {
      return sgpp::base::GridType::LinearL0Boundary;
    } else if (inputLower.compare("linearboundary") == 0) {
      return sgpp::base::GridType::LinearBoundary;
    } else if (inputLower.compare("linearstretchedboundary") == 0) {
      return sgpp::base::GridType::LinearStretchedBoundary;
    } else if (inputLower.compare("lineartruncatedboundary") == 0) {
      return sgpp::base::GridType::LinearTruncatedBoundary;
    } else if (inputLower.compare("modlinear") == 0) {
      return sgpp::base::GridType::ModLinear;
    } else if (inputLower.compare("poly") == 0) {
      return sgpp::base::GridType::Poly;
    } else if (inputLower.compare("polyboundary") == 0) {
      return sgpp::base::GridType::PolyBoundary;
    } else if (inputLower.compare("modpoly") == 0) {
      return sgpp::base::GridType::ModPoly;
    } else if (inputLower.compare("modwavelet") == 0) {
      return sgpp::base::GridType::ModWavelet;
    } else if (inputLower.compare("modbspline") == 0) {
      return sgpp::base::GridType::ModBspline;
    } else if (inputLower.compare("prewavelet") == 0) {
      return sgpp::base::GridType::Prewavelet;
    } else if (inputLower.compare("squareroot") == 0) {
      return sgpp::base::GridType::SquareRoot;
    } else if (inputLower.compare("periodic") == 0) {
      return sgpp::base::GridType::Periodic;
    } else if (inputLower.compare("linearclenshawcurtis") == 0) {
      return sgpp::base::GridType::LinearClenshawCurtis;
    } else if (inputLower.compare("bspline") == 0) {
      return sgpp::base::GridType::Bspline;
    } else if (inputLower.compare("bsplineboundary") == 0) {
      return sgpp::base::GridType::BsplineBoundary;
    } else if (inputLower.compare("bsplineclenshawcurtis") == 0) {
      return sgpp::base::GridType::BsplineClenshawCurtis;
    } else if (inputLower.compare("wavelet") == 0) {
      return sgpp::base::GridType::Wavelet;
    } else if (inputLower.compare("waveletboundary") == 0) {
      return sgpp::base::GridType::WaveletBoundary;
    } else if (inputLower.compare("fundamentalspline") == 0) {
      return sgpp::base::GridType::FundamentalSpline;
    } else if (inputLower.compare("modfundamentalspline") == 0) {
      return sgpp::base::GridType::ModFundamentalSpline;
    } else if (inputLower.compare("modbsplineclenshawcurtis") == 0) {
      return sgpp::base::GridType::ModBsplineClenshawCurtis;
    } else if (inputLower.compare("linearstencil") == 0) {
      return sgpp::base::GridType::LinearStencil;
    } else if (inputLower.compare("modlinearstencil") == 0) {
      return sgpp::base::GridType::ModLinearStencil;
    } else {
      std::string errorMsg = "Failed to convert string \"" + input + "\" to any known GridType";
      throw data_exception(errorMsg.c_str());
    }
  }
};
}  // namespace base
}  // namespace sgpp
