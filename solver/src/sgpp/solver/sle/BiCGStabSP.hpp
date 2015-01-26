// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef BICGSTABSP_HPP
#define BICGSTABSP_HPP

#include <sgpp/base/operation/OperationMatrixSP.hpp>
#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <sgpp/solver/SLESolverSP.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    class BiCGStabSP : public SLESolverSP {
      private:


      public:
        /**
         * Std-Constructor
         */
        BiCGStabSP(size_t imax, float epsilon);

        /**
         * Std-Destructor
         */
        virtual ~BiCGStabSP();

        /**
         * max_threashold is ignored in this solver
         *
         * Reference:
         * http://www.iue.tuwien.ac.at/phd/heinreichsberger/node70.html
         * http://www.numerik.math.tu-graz.ac.at/kurse/lgs/SIMNET6.pdf
         * http://netlib.org
         */
        virtual void solve(SGPP::base::OperationMatrixSP& SystemMatrix, SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& b, bool reuse = false, bool verbose = false, float max_threshold = -1.0);
    };

  }
}

#endif /* BICGSTABSP_HPP */