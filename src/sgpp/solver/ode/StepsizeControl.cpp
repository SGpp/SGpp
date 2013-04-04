/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)


#include "base/grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/StepsizeControl.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/tools/GridPrinter.hpp"
#include "base/exception/solver_exception.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
namespace sg {
  namespace solver {

    StepsizeControl::StepsizeControl(size_t imax, double timestepSize, double eps, double sc, sg::base::ScreenOutput* screen, double gamma)
      : ODESolver(imax, timestepSize), myScreen(screen), _gamma(gamma) {
      this->residuum = 0.0;
      this->myEps = eps;
      this->mySC = sc;
      useCoarsen = true; // todo

    }

    StepsizeControl::~StepsizeControl() {
    }

    double StepsizeControl::norm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2) {
      return twoNorm(System, dv1, dv2);
    }

    double StepsizeControl::twoNorm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& dv1, sg::base::DataVector& dv2) {
      dv1.sub(dv2);

      return sqrt(dv1.dotProduct(dv1));
    }

    double StepsizeControl::maxNorm(sg::pde::OperationParabolicPDESolverSystem& System, sg::base::DataVector& YkImEul, sg::base::DataVector& YkImEulOld) {
      double max = 0.0;
      double sc = this->mySC;

      double* OldData = YkImEulOld.getPointer();
      double* Data = YkImEul.getPointer();

      // calculate the max norm
      if (!useCoarsen) {
        for (size_t j = 0; j < System.getGridCoefficientsForCG()->getSize(); j++) {
          double t2 = std::max(fabs(Data[j]), fabs(OldData[j]));
          double tmpData = fabs(Data[j] - OldData[j]) / std::max(sc, t2);

          if (max < fabs(tmpData))
            max = fabs(tmpData);
        }
      } else {
        //      std::cout << "YK" << YkImEul.getSize()<< " YKold" << YkImEulOld.getSize()  << std::endl;

        sg::base::GridStorage* gs = System.getGridStorage();
        sg::base::GridStorage* ogs = System.getOldGridStorage();
        sg::base::GridStorage::grid_map_iterator q;

        for (sg::base::GridStorage::grid_map_iterator p = gs->begin(); p != gs->end(); ++p) {
          q = ogs->find(p->first);

          //   std::cout << "YK" << ((p->first)->toString() == (q->first)->toString() ) << " " << (p->first)->toString() << " YKold" << ((q->first)->equals(*p->first)) << " "<< (q->first)->toString() << std::endl;
          if ((q->first)->equals(*p->first)) {

            long unsigned int i = p->second;
            long unsigned int j = q->second;
            //  std::cout <<time<< " "<< (p->first)->toString() << " "<<i<<" " << Data[i]<< " " << (q->first)->toString() <<" "<<j<< " "<< OldData[j] << std::endl;
            double t2 = std::max(fabs(Data[i]), fabs(OldData[j]));
            double tmpData = fabs(Data[i] - OldData[j]) / std::max(sc, t2);

            if (max < fabs(tmpData)) {

              max = fabs(tmpData);
            }
          }
        }

        //std::cout << "u1 " << max << " e"<<(max >= epsilon) << std::endl;

      }

      return max;

    }
    void StepsizeControl::solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep, bool verbose) {
      size_t allIter = 0;
      sg::base::DataVector* rhs;
      sg::base::DataVector YkAdBas(System.getGridCoefficients()->getSize());
      sg::base::DataVector YkImEul(System.getGridCoefficients()->getSize());

      double tmp_timestepsize = this->myEpsilon;
      double tmp_timestepsize_old = tmp_timestepsize;
      double tmp_timestepsize_new = tmp_timestepsize;
      double epsilon = this->myEps;

      double maxTimestep = static_cast<double> (this->nMaxIterations) * tmp_timestepsize;

      size_t maxIter = this->nMaxIterations * 10000;

      double time = 0.0;

      std::ofstream fileout;

      //fileout.open(filename.c_str());
      fileout.open(filename.c_str(), std::ofstream::app); // apend to file
      fileout << std::endl;

      System.getGridCoefficientsForSC(YkImEul);

      rhs = NULL;

      for (size_t i = 0; i < maxIter && time < maxTimestep; i++) {


        YkAdBas.resize(System.getGridCoefficients()->getSize());

        YkImEul.resize(System.getGridCoefficients()->getSize());

        predictor(LinearSystemSolver, System, tmp_timestepsize, YkAdBas, YkImEul, rhs);

        corrector(LinearSystemSolver, System, tmp_timestepsize, YkImEul, rhs);

        double tmp  = norm(System, YkImEul, YkAdBas);

        tmp_timestepsize_new = nextTimestep(tmp_timestepsize, tmp_timestepsize_old, tmp, epsilon);

        if (0.8 * tmp_timestepsize > tmp_timestepsize_new) {
          if (_gamma <= 0.0)
            tmp_timestepsize = tmp_timestepsize_new;
          else
            tmp_timestepsize = _gamma * tmp_timestepsize;

          System.abortTimestep();
          allIter += LinearSystemSolver.getNumberIterations();

        } else {
          fileout << i << " " << (tmp_timestepsize_new - tmp_timestepsize) << " " << time << " " << tmp_timestepsize << std::endl;
          time += tmp_timestepsize;
          allIter += LinearSystemSolver.getNumberIterations();

          if (verbose == true) {
            if (myScreen == NULL) {
              std::cout << "Final residuum " << LinearSystemSolver.getResiduum() << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")" << std::endl;
            }
          }

          if (myScreen != NULL) {
            std::stringstream soutput;

            soutput << " Final residuum " << LinearSystemSolver.getResiduum() << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")";

            if (i < this->nMaxIterations - 1) {
              myScreen->update((size_t)(((double)(time) * 100.0) / ((double)maxTimestep)), soutput.str());
            } else {
              myScreen->update(100, soutput.str());
            }
          }

          if (bIdentifyLastStep == false) {
            System.coarsenAndRefine(false);
          } else {
            if (i < (this->nMaxIterations - 1)) {
              System.coarsenAndRefine(false);
            } else {
              System.coarsenAndRefine(true);
            }
          }

          System.saveAlpha();

          tmp_timestepsize_old = tmp_timestepsize;
          tmp_timestepsize = tmp_timestepsize_new;

          // avoid small last time steps
          if (maxTimestep - time < 1.3 * tmp_timestepsize) {
            tmp_timestepsize = maxTimestep - time;
          }

          // adapt size of last time step
          tmp_timestepsize = std::min(tmp_timestepsize, maxTimestep - time);

        }

      }

      fileout.close();

      // write some empty lines to console
      if (myScreen != NULL) {
        myScreen->writeEmptyLines(2);
      }



    }
  }


}
