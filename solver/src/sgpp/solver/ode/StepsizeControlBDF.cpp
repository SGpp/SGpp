/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/StepsizeControlBDF.hpp>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

namespace sg {
  namespace solver {

    StepsizeControlBDF::StepsizeControlBDF(size_t nTimesteps, double timestepSize, double eps, sg::base::ScreenOutput* screen)
      : VarTimestep("AdBas", "BDF2", nTimesteps, timestepSize, eps, screen) {
      std::stringstream fnsstream;
      fnsstream << "Time_" << "SCBDF" << this->myEps << ".gnuplot";
      filename = fnsstream.str();
    }

    StepsizeControlBDF::~StepsizeControlBDF() {
    }

    double StepsizeControlBDF::nextTimestep(double tmp_timestepsize, double tmp_timestepsize_old, double norm, double epsilon) {

      // double deltaY = u/(3.0*(1.0+tmp_timestepsize/tmp_timestepsize_old));
      /*
          double epsilon = 0.001;
          YkBDF2.sub(YkF23);
          double u  = std::max(YkBDF2.max(),-YkBDF2.min());
          double tD = tmp_timestepsize/tmp_timestepsize_old;
          double C1 = (1.0+1.0/tD)/6.0;
          double Cp = -(1.0+tD)*(1.0+tD)/(6*tD*(1.0+2.0*tD));
          double alpha0 = (1+2*tD)/(1+tD); //??
          double deltaY = alpha0*C1*u/(tmp_timestepsize*(C1-Cp));
          tmp_timestepsize_new = tmp_timestepsize * sqrt(epsilon/deltaY);
      */
      return tmp_timestepsize * pow(epsilon / norm, (double)1.0 / (double)3.0);

    }
  }
}
