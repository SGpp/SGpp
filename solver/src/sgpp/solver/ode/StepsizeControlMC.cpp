// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/StepsizeControlMC.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/exception/solver_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

namespace SGPP {
namespace solver {

StepsizeControlMC::StepsizeControlMC(size_t imax, float_t timestepSize, float_t eps,
                                     SGPP::base::ScreenOutput* screen)
    : VarTimestep("MPR", "CrNic", imax, timestepSize, eps, screen) {
  std::stringstream fnsstream;
  fnsstream << "Time_"
            << "SCMC" << eps << ".gnuplot";
  filename = fnsstream.str();
}

StepsizeControlMC::~StepsizeControlMC() {}

float_t StepsizeControlMC::nextTimestep(float_t tmp_timestepsize, float_t tmp_timestepsize_old,
                                        float_t norm, float_t epsilon) {
  float_t deltaY = norm;
  return tmp_timestepsize *
         std::max(0.67, std::min(1.5, pow(epsilon / deltaY, (float_t)1.0 / (float_t)3.0)));
}
}  // namespace solver
}  // namespace SGPP
