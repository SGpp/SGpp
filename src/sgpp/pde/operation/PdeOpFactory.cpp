/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "pde/operation/PdeOpFactory.hpp"

#include <cstring>

#include "exception/factory_exception.hpp"

#include "grid/type/PrewaveletGrid.hpp"

#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"
#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"
#include "basis/modlinear/operation/pde/OperationLaplaceModLinear.hpp"
#include "basis/prewavelet/operation/datadriven/OperationLaplacePrewavelet.hpp"
#include "basis/linearstretched/noboundary/operation/pde/OperationLaplaceLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/pde/OperationLaplaceLinearStretchedBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"
#include "basis/linear/boundary/operation/pde/OperationLTwoDotProductLinearBoundary.hpp"
#include "basis/linearstretched/noboundary/operation/pde/OperationLTwoDotProductLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/pde/OperationLTwoDotProductLinearStretchedBoundary.hpp"

// @todo (heinecke) remove this when done
#include "basis/linearstretched/boundary/operation/common/OperationUpDownTestLinearStretchedBoundary.hpp"

namespace sg
{

namespace GOperationFactory
{

  base::OperationMatrix* createOperationLaplace(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new pde::OperationLaplaceLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new pde::OperationLaplaceLinearBoundary(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "modlinear") == 0 )
      {
        return new pde::OperationLaplaceModLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "prewavelet") == 0 )
      {
        return new pde::OperationLaplacePrewavelet(grid.getStorage(),
                                                   ((base::PrewaveletGrid*) &grid)->getShadowStorage());
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0 )
      {
        return new pde::OperationLaplaceLinearStretched(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
      {
        return new pde::OperationLaplaceLinearStretchedBoundary(grid.getStorage());
      }
    else
      {
        throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
      }
  }

  base::OperationMatrix* createOperationLaplace(base::Grid& grid, sg::base::DataVector& coef)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new pde::OperationLaplaceLinear(grid.getStorage(), coef);
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new pde::OperationLaplaceLinearBoundary(grid.getStorage(), coef);
      }
    else
      {
        throw base::factory_exception("OperationLaplace (with coefficients) is not implemented for this grid type.");
      }
  }

  base::OperationMatrix* createOperationLTwoDotProduct(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new pde::OperationLTwoDotProductLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new pde::OperationLTwoDotProductLinearBoundary(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0)
      {
        return new pde::OperationLTwoDotProductLinearStretched(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
      {
        return new pde::OperationLTwoDotProductLinearStretchedBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationUpDownTest(base::Grid& grid)
  {
    if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
      {
        return new pde::OperationUpDownTestLinearStretchedBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
  }


}
}

