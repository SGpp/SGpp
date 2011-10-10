/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "datadriven/operation/DatadrivenOpFactory.hpp"

#include <cstring>

#include "exception/factory_exception.hpp"

#include "grid/type/PolyGrid.hpp"
#include "grid/type/ModPolyGrid.hpp"
#include "grid/type/PrewaveletGrid.hpp"
#include "grid/type/ModBsplineGrid.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationTestLinear.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"
#include "basis/modbspline/operation/datadriven/OperationTestModBspline.hpp"
#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"
#include "basis/poly/operation/datadriven/OperationTestPoly.hpp"
#include "basis/modpoly/operation/datadriven/OperationTestModPoly.hpp"
#include "basis/modwavelet/operation/datadriven/OperationTestModWavelet.hpp"
#include "basis/prewavelet/operation/datadriven/OperationTestPrewavelet.hpp"
#include "basis/linearstretched/boundary/operation/datadriven/OperationTestLinearStretchedBoundary.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationTestLinearStretched.hpp"


namespace sg
{

namespace op_factory
{

  datadriven::OperationTest* createOperationTest(base::Grid& grid)
  {
    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new datadriven::OperationTestLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new datadriven::OperationTestLinearBoundary(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "modBspline") == 0 )
      {
        return new datadriven::OperationTestModBspline(grid.getStorage(),
                                                   ((base::ModBsplineGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modlinear") == 0 )
      {
        return new datadriven::OperationTestModLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "poly") == 0 )
      {
        return new datadriven::OperationTestPoly(grid.getStorage(),
                                             ((base::PolyGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modpoly") == 0 )
      {
        return new datadriven::OperationTestModPoly(grid.getStorage(),
                                                ((base::ModPolyGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modWavelet") == 0 )
      {
        return new datadriven::OperationTestModWavelet(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "prewavelet") == 0 )
      {
        return new datadriven::OperationTestPrewavelet(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0 )
      {
        return new datadriven::OperationTestLinearStretched(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
      {
        return new datadriven::OperationTestLinearStretchedBoundary(grid.getStorage());
      }

    else
      throw base::factory_exception("OperationTest is not implemented for this grid type.");
  }


}
}

