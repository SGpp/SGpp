// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/StepsizeControlBDF.hpp>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    StepsizeControlBDF::StepsizeControlBDF(size_t nTimesteps, float_t timestepSize, float_t eps, SGPP::base::ScreenOutput* screen)
      : VarTimestep("AdBas", "BDF2", nTimesteps, timestepSize, eps, screen) {
      std::stringstream fnsstream;
      fnsstream << "Time_" << "SCBDF" << this->myEps << ".gnuplot";
      filename = fnsstream.str();
    }

    StepsizeControlBDF::~StepsizeControlBDF() {
    }

    float_t StepsizeControlBDF::nextTimestep(float_t tmp_timestepsize, float_t tmp_timestepsize_old, float_t norm, float_t epsilon) {

      // float_t deltaY = u/(3.0*(1.0+tmp_timestepsize/tmp_timestepsize_old));
      /*
          float_t epsilon = 0.001;
          YkBDF2.sub(YkF23);
          float_t u  = std::max(YkBDF2.max(),-YkBDF2.min());
          float_t tD = tmp_timestepsize/tmp_timestepsize_old;
          float_t C1 = (1.0+1.0/tD)/6.0;
          float_t Cp = -(1.0+tD)*(1.0+tD)/(6*tD*(1.0+2.0*tD));
          float_t alpha0 = (1+2*tD)/(1+tD); //??
          float_t deltaY = alpha0*C1*u/(tmp_timestepsize*(C1-Cp));
          tmp_timestepsize_new = tmp_timestepsize * sqrt(epsilon/deltaY);
      */
      return tmp_timestepsize * pow(epsilon / norm, (float_t)1.0 / (float_t)3.0);

    }
  }
}