/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BICGSTABSP_HPP
#define BICGSTABSP_HPP

#include "base/operation/OperationMatrixSP.hpp"
#include "base/datatypes/DataVectorSP.hpp"

#include "solver/SLESolverSP.hpp"

#include <iostream>

namespace sg {
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
        virtual void solve(sg::base::OperationMatrixSP& SystemMatrix, sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& b, bool reuse = false, bool verbose = false, float max_threshold = -1.0);
    };

  }
}

#endif /* BICGSTABSP_HPP */
