/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/operation/BaseOpFactory.hpp"

#include "base/exception/factory_exception.hpp"

#include "base/grid/type/PolyGrid.hpp"
#include "base/grid/type/ModPolyGrid.hpp"
#include "base/grid/type/PrewaveletGrid.hpp"
#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"
#include "base/grid/type/BsplineGrid.hpp"
#include "base/grid/type/BsplineTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModBsplineGrid.hpp"
#include "base/grid/type/LinearClenshawCurtisGrid.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"
#include "base/grid/type/WaveletGrid.hpp"
#include "base/grid/type/WaveletTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModWaveletGrid.hpp"

#include "base/basis/linear/noboundary/operation/OperationStencilHierarchisationLinear.hpp"
#include "base/basis/linear/modified/operation/OperationStencilHierarchisationModLinear.hpp"
#include "base/basis/linear/noboundary/operation/OperationHierarchisationLinear.hpp"
#include "base/basis/linear/modified/operation/OperationHierarchisationModLinear.hpp"
#include "base/basis/linear/boundary/operation/OperationHierarchisationLinearBoundary.hpp"
#include "base/basis/linearstretched/noboundary/operation/OperationHierarchisationLinearStretched.hpp"
#include "base/basis/linearstretched/boundary/operation/OperationHierarchisationLinearStretchedBoundary.hpp"
#include "base/basis/poly/operation/OperationHierarchisationPoly.hpp"
#include "base/basis/modpoly/operation/OperationHierarchisationModPoly.hpp"
#include "base/basis/prewavelet/operation/OperationHierarchisationPrewavelet.hpp"

#include "base/operation/OperationQuadrature.hpp"
#include "base/basis/linear/noboundary/operation/OperationQuadratureLinear.hpp"
#include "base/basis/linear/noboundary/operation/OperationFirstMomentLinear.hpp"
#include "base/basis/linear/noboundary/operation/OperationSecondMomentLinear.hpp"
#include "base/basis/linear/boundary/operation/OperationQuadratureLinearBoundary.hpp"
#include "base/basis/poly/operation/OperationQuadraturePoly.hpp"

#include "base/basis/prewavelet/operation/OperationConvertPrewavelet.hpp"

#include "base/basis/linear/noboundary/operation/OperationEvalLinear.hpp"
#include "base/basis/linear/boundary/operation/OperationEvalLinearBoundary.hpp"
#include "base/basis/linear/modified/operation/OperationEvalModLinear.hpp"
#include "base/basis/poly/operation/OperationEvalPoly.hpp"
#include "base/basis/modpoly/operation/OperationEvalModPoly.hpp"
#include "base/basis/bspline/modified/operation/OperationEvalModBspline.hpp"
#include "base/basis/wavelet/modified/operation/OperationEvalModWavelet.hpp"
#include "base/basis/prewavelet/operation/OperationEvalPrewavelet.hpp"
#include "base/basis/linearstretched/noboundary/operation/OperationEvalLinearStretched.hpp"
#include "base/basis/linearstretched/boundary/operation/OperationEvalLinearStretchedBoundary.hpp"

#include "base/basis/linear/noboundary/operation/OperationMultipleEvalLinear.hpp"
#include "base/basis/linear/boundary/operation/OperationMultipleEvalLinearBoundary.hpp"
#include "base/basis/linear/modified/operation/OperationMultipleEvalModLinear.hpp"
#include "base/basis/poly/operation/OperationMultipleEvalPoly.hpp"
#include "base/basis/modpoly/operation/OperationMultipleEvalModPoly.hpp"
#include "base/basis/bspline/modified/operation/OperationMultipleEvalModBspline.hpp"
#include "base/basis/wavelet/modified/operation/OperationMultipleEvalModWavelet.hpp"
#include "base/basis/prewavelet/operation/OperationMultipleEvalPrewavelet.hpp"
#include "base/basis/linearstretched/noboundary/operation/OperationMultipleEvalLinearStretched.hpp"
#include "base/basis/linearstretched/boundary/operation/OperationMultipleEvalLinearStretchedBoundary.hpp"

#include "base/basis/bspline/boundary/operation/OperationNaiveEvalBsplineBoundary.hpp"
#include "base/basis/bspline/clenshawcurtis/operation/OperationNaiveEvalBsplineClenshawCurtis.hpp"
#include "base/basis/bspline/modified/operation/OperationNaiveEvalModBspline.hpp"
#include "base/basis/bspline/noboundary/operation/OperationNaiveEvalBspline.hpp"
#include "base/basis/linear/boundary/operation/OperationNaiveEvalLinearBoundary.hpp"
#include "base/basis/linear/clenshawcurtis/operation/OperationNaiveEvalLinearClenshawCurtis.hpp"
#include "base/basis/linear/modified/operation/OperationNaiveEvalModLinear.hpp"
#include "base/basis/linear/noboundary/operation/OperationNaiveEvalLinear.hpp"
#include "base/basis/wavelet/boundary/operation/OperationNaiveEvalWaveletBoundary.hpp"
#include "base/basis/wavelet/modified/operation/OperationNaiveEvalModWavelet.hpp"
#include "base/basis/wavelet/noboundary/operation/OperationNaiveEvalWavelet.hpp"

#include "base/basis/bspline/boundary/operation/OperationNaiveEvalGradientBsplineBoundary.hpp"
#include "base/basis/bspline/clenshawcurtis/operation/OperationNaiveEvalGradientBsplineClenshawCurtis.hpp"
#include "base/basis/bspline/modified/operation/OperationNaiveEvalGradientModBspline.hpp"
#include "base/basis/bspline/noboundary/operation/OperationNaiveEvalGradientBspline.hpp"
#include "base/basis/wavelet/boundary/operation/OperationNaiveEvalGradientWaveletBoundary.hpp"
#include "base/basis/wavelet/modified/operation/OperationNaiveEvalGradientModWavelet.hpp"
#include "base/basis/wavelet/noboundary/operation/OperationNaiveEvalGradientWavelet.hpp"

#include "base/basis/bspline/boundary/operation/OperationNaiveEvalHessianBsplineBoundary.hpp"
#include "base/basis/bspline/clenshawcurtis/operation/OperationNaiveEvalHessianBsplineClenshawCurtis.hpp"
#include "base/basis/bspline/modified/operation/OperationNaiveEvalHessianModBspline.hpp"
#include "base/basis/bspline/noboundary/operation/OperationNaiveEvalHessianBspline.hpp"
#include "base/basis/wavelet/boundary/operation/OperationNaiveEvalHessianWaveletBoundary.hpp"
#include "base/basis/wavelet/modified/operation/OperationNaiveEvalHessianModWavelet.hpp"
#include "base/basis/wavelet/noboundary/operation/OperationNaiveEvalHessianWavelet.hpp"

#include "base/basis/bspline/boundary/operation/OperationNaiveEvalPartialDerivativeBsplineBoundary.hpp"
#include "base/basis/bspline/clenshawcurtis/operation/OperationNaiveEvalPartialDerivativeBsplineClenshawCurtis.hpp"
#include "base/basis/bspline/modified/operation/OperationNaiveEvalPartialDerivativeModBspline.hpp"
#include "base/basis/bspline/noboundary/operation/OperationNaiveEvalPartialDerivativeBspline.hpp"
#include "base/basis/wavelet/boundary/operation/OperationNaiveEvalPartialDerivativeWaveletBoundary.hpp"
#include "base/basis/wavelet/modified/operation/OperationNaiveEvalPartialDerivativeModWavelet.hpp"
#include "base/basis/wavelet/noboundary/operation/OperationNaiveEvalPartialDerivativeWavelet.hpp"

#include <cstring>

namespace sg {

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
      } else
        throw base::factory_exception("OperationEval is not implemented for this grid type.");
    }

    base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid, base::DataMatrix* dataset) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new base::OperationMultipleEvalLinear(grid.getStorage(), dataset);
      } else if (strcmp(grid.getType(), "linearBoundary") == 0
                 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new base::OperationMultipleEvalLinearBoundary(grid.getStorage(), dataset);
      } else if (strcmp(grid.getType(), "modlinear") == 0 ) {
        return new base::OperationMultipleEvalModLinear(grid.getStorage(), dataset);
      } else if (strcmp(grid.getType(), "poly") == 0 ) {
        return new base::OperationMultipleEvalPoly(grid.getStorage(),
               ((base::PolyGrid*) &grid)->getDegree(), dataset);
      } else if (strcmp(grid.getType(), "modpoly") == 0 ) {
        return new base::OperationMultipleEvalModPoly(grid.getStorage(),
               ((base::ModPolyGrid*) &grid)->getDegree(), dataset);
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationMultipleEvalModBspline(grid.getStorage(),
               ((base::ModBsplineGrid*) &grid)->getDegree(), dataset);
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationMultipleEvalModWavelet(grid.getStorage(), dataset);
      } else if (strcmp(grid.getType(), "prewavelet") == 0 ) {
        return new base::OperationMultipleEvalPrewavelet(grid.getStorage(), dataset);
      } else if (strcmp(grid.getType(), "linearStretched") == 0 ) {
        return new base::OperationMultipleEvalLinearStretched(grid.getStorage(), dataset);
      } else if (strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 ) {
        return new base::OperationMultipleEvalLinearStretchedBoundary(grid.getStorage(), dataset);
      }

      else
        throw base::factory_exception("OperationMultipleEval is not implemented for this grid type.");
    }

    base::OperationNaiveEval* createOperationNaiveEval(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0 ) {
        return new base::OperationNaiveEvalLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "modlinear") == 0 ) {
        return new base::OperationNaiveEvalModLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalLinearBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearClenshawCurtis") == 0 ) {
        return new base::OperationNaiveEvalLinearClenshawCurtis(grid.getStorage(),
               dynamic_cast<base::LinearClenshawCurtisGrid &>(grid).getCosineTable());
      } else if (strcmp(grid.getType(), "Bspline") == 0 ) {
        return new base::OperationNaiveEvalBspline(grid.getStorage(),
               dynamic_cast<base::BsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationNaiveEvalModBspline(grid.getStorage(),
               dynamic_cast<base::ModBsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalBsplineBoundary(grid.getStorage(),
               dynamic_cast<base::BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0 ) {
        return new base::OperationNaiveEvalBsplineClenshawCurtis(grid.getStorage(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getDegree(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getCosineTable());
      } else if (strcmp(grid.getType(), "Wavelet") == 0 ) {
        return new base::OperationNaiveEvalWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationNaiveEvalModWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalWaveletBoundary(grid.getStorage());
      } else
        throw base::factory_exception("OperationNaiveEval is not implemented for this grid type.");
    }

    base::OperationNaiveEvalGradient* createOperationNaiveEvalGradient(base::Grid& grid) {

      if (strcmp(grid.getType(), "Bspline") == 0 ) {
        return new base::OperationNaiveEvalGradientBspline(grid.getStorage(),
               dynamic_cast<base::BsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationNaiveEvalGradientModBspline(grid.getStorage(),
               dynamic_cast<base::ModBsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalGradientBsplineBoundary(grid.getStorage(),
               dynamic_cast<base::BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0 ) {
        return new base::OperationNaiveEvalGradientBsplineClenshawCurtis(grid.getStorage(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getDegree(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getCosineTable());
      } else if (strcmp(grid.getType(), "Wavelet") == 0 ) {
        return new base::OperationNaiveEvalGradientWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationNaiveEvalGradientModWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalGradientWaveletBoundary(grid.getStorage());
      } else
        throw base::factory_exception("OperationNaiveEvalGradient is not implemented for this grid type.");
    }

    base::OperationNaiveEvalHessian* createOperationNaiveEvalHessian(base::Grid& grid) {

      if (strcmp(grid.getType(), "Bspline") == 0 ) {
        return new base::OperationNaiveEvalHessianBspline(grid.getStorage(),
               dynamic_cast<base::BsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationNaiveEvalHessianModBspline(grid.getStorage(),
               dynamic_cast<base::ModBsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalHessianBsplineBoundary(grid.getStorage(),
               dynamic_cast<base::BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0 ) {
        return new base::OperationNaiveEvalHessianBsplineClenshawCurtis(grid.getStorage(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getDegree(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getCosineTable());
      } else if (strcmp(grid.getType(), "Wavelet") == 0 ) {
        return new base::OperationNaiveEvalHessianWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationNaiveEvalHessianModWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalHessianWaveletBoundary(grid.getStorage());
      } else
        throw base::factory_exception("OperationNaiveEvalHessian is not implemented for this grid type.");
    }

    base::OperationNaiveEvalPartialDerivative* createOperationNaiveEvalPartialDerivative(base::Grid& grid) {

      if (strcmp(grid.getType(), "Bspline") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeBspline(grid.getStorage(),
               dynamic_cast<base::BsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "modBspline") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeModBspline(grid.getStorage(),
               dynamic_cast<base::ModBsplineGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeBsplineBoundary(grid.getStorage(),
               dynamic_cast<base::BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
      } else if (strcmp(grid.getType(), "BsplineClenshawCurtis") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeBsplineClenshawCurtis(grid.getStorage(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getDegree(),
               dynamic_cast<base::BsplineClenshawCurtisGrid &>(grid).getCosineTable());
      } else if (strcmp(grid.getType(), "Wavelet") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "modWavelet") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeModWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0 ) {
        return new base::OperationNaiveEvalPartialDerivativeWaveletBoundary(grid.getStorage());
      } else
        throw base::factory_exception("OperationNaiveEvalPartialDerivative is not implemented for this grid type.");
    }

  }
}

