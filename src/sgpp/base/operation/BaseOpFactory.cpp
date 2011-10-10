/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/operation/BaseOpFactory.hpp"

#include <cstring>

#include "exception/factory_exception.hpp"

#include "grid/type/PolyGrid.hpp"
#include "grid/type/ModPolyGrid.hpp"
#include "grid/type/PrewaveletGrid.hpp"
#include "grid/type/ModBsplineGrid.hpp"

#include "basis/linear/noboundary/operation/common/OperationHierarchisationLinear.hpp"
#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationHierarchisationLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/common/OperationHierarchisationLinearStretchedBoundary.hpp"
#include "basis/poly/operation/common/OperationHierarchisationPoly.hpp"
#include "basis/modpoly/operation/common/OperationHierarchisationModPoly.hpp"
#include "basis/prewavelet/operation/common/OperationHierarchisationPrewavelet.hpp"
#include "basis/modbspline/operation/common/OperationHierarchisationModBspline.hpp"
#include "basis/modwavelet/operation/common/OperationHierarchisationModWavelet.hpp"

#include "operation/common/OperationQuadrature.hpp"
#include "basis/linear/noboundary/operation/common/OperationQuadratureLinear.hpp"
#include "basis/poly/operation/common/OperationQuadraturePoly.hpp"

#include "basis/prewavelet/operation/common/OperationConvertPrewavelet.hpp"

#include "basis/linear/noboundary/operation/common/OperationEvalLinear.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/modlinear/operation/common/OperationEvalModLinear.hpp"
#include "basis/poly/operation/common/OperationEvalPoly.hpp"
#include "basis/modpoly/operation/common/OperationEvalModPoly.hpp"
#include "basis/modbspline/operation/common/OperationEvalModBspline.hpp"
#include "basis/modwavelet/operation/common/OperationEvalModWavelet.hpp"
#include "basis/prewavelet/operation/common/OperationEvalPrewavelet.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationEvalLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/common/OperationEvalLinearStretchedBoundary.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalLinear.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationMultipleEvalLinearBoundary.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalModLinear.hpp"
#include "basis/poly/operation/datadriven/OperationMultipleEvalPoly.hpp"
#include "basis/modpoly/operation/datadriven/OperationMultipleEvalModPoly.hpp"
#include "basis/modbspline/operation/datadriven/OperationMultipleEvalModBspline.hpp"
#include "basis/modwavelet/operation/datadriven/OperationMultipleEvalModWavelet.hpp"
#include "basis/prewavelet/operation/datadriven/OperationMultipleEvalPrewavelet.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalLinearStretched.hpp"
#include "basis/linearstretched/boundary/operation/datadriven/OperationMultipleEvalLinearStretchedBoundary.hpp"


namespace sg
{

namespace op_factory
{

  base::OperationHierarchisation* createOperationHierarchisation(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new base::OperationHierarchisationLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "modlinear") == 0 )
      {
        return new base::OperationHierarchisationModLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
            || strcmp(grid.getType(), "TruncatedTrapezoid") == 0
            || strcmp(grid.getType(), "squareRoot") == 0)
      {
        return new base::OperationHierarchisationLinearBoundary(grid.getStorage());
      }

    else if(strcmp(grid.getType(), "linearStretched") == 0 )
      {
        return new base::OperationHierarchisationLinearStretched(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
      {
        return new base::OperationHierarchisationLinearStretchedBoundary(grid.getStorage());
      }
//    else if(strcmp(grid.getType(), "poly") == 0 )
//      {
//        return new base::OperationHierarchisationPoly(grid.getStorage(),
//                                                      ((base::PolyGrid*) &grid)->getDegree());
//      }
    else if(strcmp(grid.getType(), "modpoly") == 0 )
      {
        return new base::OperationHierarchisationModPoly(grid.getStorage(),
                                                         ((base::ModPolyGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "prewavelet") == 0 )
      {
        return new base::OperationHierarchisationPrewavelet(grid.getStorage(),
                                                      ((base::PrewaveletGrid*) &grid)->getShadowStorage());
      }
    else if(strcmp(grid.getType(), "modBspline") == 0 )
      {
        return new base::OperationHierarchisationModBspline(grid.getStorage(),
                                                      ((base::ModBsplineGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modWavelet") == 0 )
      {
        return new base::OperationHierarchisationModWavelet(grid.getStorage());
      }

    else
      throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
  }

  base::OperationQuadrature* createOperationQuadrature(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new base::OperationQuadratureLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "poly") == 0 )
      {
        if(((base::PolyGrid*) &grid)->getDegree()>3) {
          throw base::factory_exception("OperationQuadrature is not implemented for polynomials with degree higher than 3.");
        }
        else {
          return new base::OperationQuadraturePoly(grid.getStorage(), ((base::PolyGrid*) &grid)->getDegree());
        }
      }
    else
      throw base::factory_exception("OperationQuadrature is not implemented for this grid type.");
  }

  base::OperationConvert* createOperationConvert(base::Grid& grid)
  {
    if(strcmp(grid.getType(), "prewavelet") == 0 )
      {
        return new base::OperationConvertPrewavelet(grid.getStorage(),
                                                    ((base::PrewaveletGrid*) &grid)->getShadowStorage());
      }

    else
      throw base::factory_exception("OperationConvert is not implemented for this grid type.");
  }

  base::OperationMatrix* createOperationIdentity(base::Grid& grid)
  {
    return new base::OperationIdentity();
  }

  base::OperationEval* createOperationEval(base::Grid& grid)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new base::OperationEvalLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
            || strcmp(grid.getType(), "TruncatedTrapezoid") == 0
            || strcmp(grid.getType(), "squareRoot") == 0)
      {
        return new base::OperationEvalLinearBoundary(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "modlinear") == 0 )
      {
        return new base::OperationEvalModLinear(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "poly") == 0 )
      {
        return new base::OperationEvalPoly(grid.getStorage(),
                                     ((base::PolyGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modpoly") == 0 )
      {
        return new base::OperationEvalModPoly(grid.getStorage(),
                                        ((base::ModPolyGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modBspline") == 0 )
      {
        return new base::OperationEvalModBspline(grid.getStorage(),
                                           ((base::ModBsplineGrid*) &grid)->getDegree());
      }
    else if(strcmp(grid.getType(), "modWavelet") == 0 )
      {
        return new base::OperationEvalModWavelet(grid.getStorage());
      }

    else if(strcmp(grid.getType(), "prewavelet") == 0 )
      {
        return new base::OperationEvalPrewavelet(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0 )
      {
        return new base::OperationEvalLinearStretched(grid.getStorage());
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
      {
        return new base::OperationEvalLinearStretchedBoundary(grid.getStorage());
      }
    else
      throw base::factory_exception("OperationEval is not implemented for this grid type.");
  }

  base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid, base::DataMatrix* dataset)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        return new base::OperationMultipleEvalLinear(grid.getStorage(), dataset);
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        return new base::OperationMultipleEvalLinearBoundary(grid.getStorage(), dataset);
      }
    else if(strcmp(grid.getType(), "modlinear") == 0 )
      {
        return new base::OperationMultipleEvalModLinear(grid.getStorage(), dataset);
      }
    else if(strcmp(grid.getType(), "poly") == 0 )
      {
        return new base::OperationMultipleEvalPoly(grid.getStorage(),
                                             ((base::PolyGrid*) &grid)->getDegree(), dataset);
      }
    else if(strcmp(grid.getType(), "modpoly") == 0 )
      {
        return new base::OperationMultipleEvalModPoly(grid.getStorage(),
                                                ((base::ModPolyGrid*) &grid)->getDegree(), dataset);
      }
    else if(strcmp(grid.getType(), "modBspline") == 0 )
      {
        return new base::OperationMultipleEvalModBspline(grid.getStorage(),
                                                   ((base::ModBsplineGrid*) &grid)->getDegree(), dataset);
      }
    else if(strcmp(grid.getType(), "modWavelet") == 0 )
      {
        return new base::OperationMultipleEvalModWavelet(grid.getStorage(), dataset);
      }
    else if(strcmp(grid.getType(), "prewavelet") == 0 )
      {
        return new base::OperationMultipleEvalPrewavelet(grid.getStorage(), dataset);
      }
    else if(strcmp(grid.getType(), "linearStretched") == 0 )
      {
        return new base::OperationMultipleEvalLinearStretched(grid.getStorage(), dataset);
      }
    else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
      {
        return new base::OperationMultipleEvalLinearStretchedBoundary(grid.getStorage(), dataset);
      }

    else
      throw base::factory_exception("OperationMultipleEval is not implemented for this grid type.");
  }



}
}

