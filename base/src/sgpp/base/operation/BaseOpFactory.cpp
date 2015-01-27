// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>

#include <sgpp/base/operation/hash/OperationStencilHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/OperationStencilHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPoly.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModPoly.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModBspline.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModWavelet.hpp>

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureLinear.hpp>
#include <sgpp/base/operation/hash/OperationFirstMomentLinear.hpp>
#include <sgpp/base/operation/hash/OperationSecondMomentLinear.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationQuadraturePoly.hpp>

#include <sgpp/base/operation/hash/OperationConvertPrewavelet.hpp>

#include <sgpp/base/operation/hash/OperationEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationEvalModPoly.hpp>
#include <sgpp/base/operation/hash/OperationEvalModBspline.hpp>
#include <sgpp/base/operation/hash/OperationEvalModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationEvalPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalPeriodic.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModPoly.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModBspline.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp>

#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {

    base::OperationHierarchisation* createOperationHierarchisation(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationHierarchisationLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearstencil") == 0 ) {
        return new base::OperationStencilHierarchisationLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "modlinearstencil") == 0 ) {
        return new base::OperationStencilHierarchisationModLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "modlinear") == 0 ) {
        return new base::OperationHierarchisationModLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0
                 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
                 || strcmp(grid.getType(), "TruncatedTrapezoid") == 0
                 || strcmp(grid.getType(), "squareRoot") == 0) {
        return new base::OperationHierarchisationLinearBoundary(grid.getStorage());
      }

      else if (strcmp(grid.getType(), "linearStretched") == 0 ) {
        return new base::OperationHierarchisationLinearStretched(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 ) {
        return new base::OperationHierarchisationLinearStretchedBoundary(grid.getStorage());
      }
      //    else if(strcmp(grid.getType(), "poly") == 0 )
      //      {
      //        return new base::OperationHierarchisationPoly(grid.getStorage(),
      //                                                      ((base::PolyGrid*) &grid)->getDegree());
      //      }
      else if (strcmp(grid.getType(), "modpoly") == 0 ) {
        return new base::OperationHierarchisationModPoly(grid.getStorage(),
               ((base::ModPolyGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "prewavelet") == 0 ) {
        return new base::OperationHierarchisationPrewavelet(grid.getStorage(),
               ((base::PrewaveletGrid*) &grid)->getShadowStorage());
      }
      //    else if(strcmp(grid.getType(), "modBspline") == 0 )
      //      {
      //        return new base::OperationHierarchisationModBspline(grid.getStorage(),
      //                                                      ((base::ModBsplineGrid*) &grid)->getDegree());
      //      }
      //    else if(strcmp(grid.getType(), "modWavelet") == 0 )
      //      {
      //        return new base::OperationHierarchisationModWavelet(grid.getStorage());
      //      }

      else
        throw base::factory_exception("OperationHierarchisation is not implemented for this grid type.");
    }

    base::OperationQuadrature* createOperationQuadrature(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationQuadratureLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0 ||
                 strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new base::OperationQuadratureLinearBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "poly") == 0 ) {
        if (((base::PolyGrid*) &grid)->getDegree() > 3) {
          throw base::factory_exception("OperationQuadrature is not implemented for polynomials with degree higher than 3.");
        } else {
          return new base::OperationQuadraturePoly(grid.getStorage(), ((base::PolyGrid*) &grid)->getDegree());
        }
      } else
        throw base::factory_exception("OperationQuadrature is not implemented for this grid type.");
    }

    base::OperationFirstMoment* createOperationFirstMoment(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationFirstMomentLinear(grid.getStorage());
      } else
        throw base::factory_exception("OperationFirstMoment is not implemented for this grid type.");
    }

    base::OperationSecondMoment* createOperationSecondMoment(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationSecondMomentLinear(grid.getStorage());
      } else
        throw base::factory_exception("OperationSecondMoment is not implemented for this grid type.");
    }

    base::OperationConvert* createOperationConvert(base::Grid& grid) {
      if (strcmp(grid.getType(), "prewavelet") == 0 ) {
        return new base::OperationConvertPrewavelet(grid.getStorage(),
               ((base::PrewaveletGrid*) &grid)->getShadowStorage());
      }

      else
        throw base::factory_exception("OperationConvert is not implemented for this grid type.");
    }

    base::OperationMatrix* createOperationIdentity(base::Grid& grid) {
      return new base::OperationIdentity();
    }

    base::OperationEval* createOperationEval(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationEvalLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0
                 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
                 || strcmp(grid.getType(), "TruncatedTrapezoid") == 0
                 || strcmp(grid.getType(), "squareRoot") == 0) {
        return new base::OperationEvalLinearBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "modlinear") == 0 ) {
        return new base::OperationEvalModLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "poly") == 0 ) {
        return new base::OperationEvalPoly(grid.getStorage(),
                                           ((base::PolyGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "modpoly") == 0 ) {
        return new base::OperationEvalModPoly(grid.getStorage(),
                                              ((base::ModPolyGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationEvalModBspline(grid.getStorage(),
               ((base::ModBsplineGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationEvalModWavelet(grid.getStorage());
      }

      else if (strcmp(grid.getType(), "prewavelet") == 0 ) {
        return new base::OperationEvalPrewavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearStretched") == 0 ) {
        return new base::OperationEvalLinearStretched(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 ) {
        return new base::OperationEvalLinearStretchedBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(),"periodic") == 0) {
    	return new base::OperationEvalPeriodic(grid.getStorage());
      } else
        throw base::factory_exception("OperationEval is not implemented for this grid type.");
    }

    base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid, base::DataMatrix &dataset) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationMultipleEvalLinear(grid, dataset);
      } else if (strcmp(grid.getType(), "linearBoundary") == 0
                 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new base::OperationMultipleEvalLinearBoundary(grid, dataset);
      } else if (strcmp(grid.getType(), "modlinear") == 0 ) {
        return new base::OperationMultipleEvalModLinear(grid, dataset);
      } else if (strcmp(grid.getType(), "poly") == 0 ) {
        return new base::OperationMultipleEvalPoly(grid,
               ((base::PolyGrid*) &grid)->getDegree(), dataset);
      } else if (strcmp(grid.getType(), "modpoly") == 0 ) {
        return new base::OperationMultipleEvalModPoly(grid,
               ((base::ModPolyGrid*) &grid)->getDegree(), dataset);
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationMultipleEvalModBspline(grid,
               ((base::ModBsplineGrid*) &grid)->getDegree(), dataset);
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationMultipleEvalModWavelet(grid, dataset);
      } else if (strcmp(grid.getType(), "prewavelet") == 0 ) {
        return new base::OperationMultipleEvalPrewavelet(grid, dataset);
      } else if (strcmp(grid.getType(), "linearStretched") == 0 ) {
        return new base::OperationMultipleEvalLinearStretched(grid, dataset);
      } else if (strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 ) {
        return new base::OperationMultipleEvalLinearStretchedBoundary(grid, dataset);
      }else if (strcmp(grid.getType(),"periodic") == 0) {
      	return new base::OperationMultipleEvalPeriodic(grid, dataset);
      }

      else
        throw base::factory_exception("OperationMultipleEval is not implemented for this grid type.");
    }



  }
}
