/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "finance/operation/FinanceOpFactory.hpp"

#include <cstring>

#include "exception/factory_exception.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLinear.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLogLinearStretchedBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLogLinear.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLogLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLinearStretchedBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLogLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLBLinear.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLBLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLELinear.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLELinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLDLinear.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLDLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLFLinear.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLFLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLinear.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLinearStretchedBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLogLinear.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLogLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLogLinearStretchedBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLogLinearBoundary.hpp"


namespace sg
{

namespace op_factory
{

  base::OperationMatrix* createOperationGamma(base::Grid& grid, base::DataMatrix& coef)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationGammaLinear(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0)
      {
        return new finance::OperationGammaLinearStretched(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
      {
        return new finance::OperationGammaLinearStretchedBoundary(grid.getStorage(), coef);
      }

    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationGammaLinearBoundary(grid.getStorage(), coef);
      }
    else
      throw base::factory_exception("OperationGamma is not implemented for this grid type.");

  }

  base::OperationMatrix* createOperationGammaLog(base::Grid& grid, base::DataMatrix& coef)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationGammaLogLinear(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0)
      {
        return new finance::OperationGammaLogLinearStretched(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
      {
        return new finance::OperationGammaLogLinearStretchedBoundary(grid.getStorage(), coef);
      }

    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationGammaLogLinearBoundary(grid.getStorage(), coef);
      }
    else
      throw base::factory_exception("OperationGammaLog is not implemented for this grid type.");

  }


  base::OperationMatrix* createOperationLB(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationLBLinear(grid.getStorage());
      }


    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationLBLinearBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationLB is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationLE(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationLELinear(grid.getStorage());
      }


    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationLELinearBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationLE is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationLD(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationLDLinear(grid.getStorage());
      }


    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationLDLinearBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationLD is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationLF(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationLFLinear(grid.getStorage());
      }


    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationLFLinearBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationLF is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationDelta(base::Grid& grid, base::DataVector& coef)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationDeltaLinear(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0)
      {
        return new finance::OperationDeltaLinearStretched(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
      {
        return new finance::OperationDeltaLinearStretchedBoundary(grid.getStorage(), coef);
      }

    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new finance::OperationDeltaLinearBoundary(grid.getStorage(), coef);
      }

    else
      throw base::factory_exception("OperationDelta is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationDeltaLog(base::Grid& grid, base::DataVector& coef)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new finance::OperationDeltaLogLinear(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0)
      {
        return new finance::OperationDeltaLogLinearStretched(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
      {
        return new finance::OperationDeltaLogLinearStretchedBoundary(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0 )
      {
        return new finance::OperationDeltaLogLinearBoundary(grid.getStorage(), coef);
      }
    else
      throw base::factory_exception("OperationDeltaLog is not implemented for this grid type.");
  }


}
}

