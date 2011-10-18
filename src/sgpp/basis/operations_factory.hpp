/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONS_FACTORY_HPP
#define OPERATIONS_FACTORY_HPP

#include <cstring>

#include "exception/factory_exception.hpp"

#include "grid/Grid.hpp"
#include "grid/type/PrewaveletGrid.hpp"
#include "grid/type/PolyGrid.hpp"
#include "grid/type/ModPolyGrid.hpp"
#include "grid/type/ModBsplineGrid.hpp"

#include "operation/common/OperationMatrix.hpp"
#include "operation/common/OperationEval.hpp"
#include "operation/datadriven/OperationMultipleEval.hpp"
#include "operation/pde/OperationParabolicPDESolverSystem.hpp"
#include "operation/datadriven/OperationTest.hpp"
#include "operation/datadriven/OperationMultipleEvalVectorized.hpp"
#include "operation/datadriven/OperationMultipleEvalVectorizedSP.hpp"
#include "operation/common/OperationHierarchisation.hpp"
#include "operation/common/OperationQuadrature.hpp"
#include "operation/common/OperationConvert.hpp"
#include "operation/common/OperationIdentity.hpp"

#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"
#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"
#include "basis/modlinear/operation/pde/OperationLaplaceModLinear.hpp"
#include "basis/prewavelet/operation/datadriven/OperationLaplacePrewavelet.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLinear.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLinearBoundary.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLogLinear.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLogLinearBoundary.hpp"

#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLELinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLBLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLFLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLDLinearBoundary.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLELinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLBLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLFLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLDLinear.hpp"

#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"
#include "basis/linear/boundary/operation/pde/OperationLTwoDotProductLinearBoundary.hpp"

#include "basis/modwavelet/operation/common/OperationEvalModWavelet.hpp"
#include "basis/linear/noboundary/operation/common/OperationEvalLinear.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/modbspline/operation/common/OperationEvalModBspline.hpp"
#include "basis/modlinear/operation/common/OperationEvalModLinear.hpp"
#include "basis/modpoly/operation/common/OperationEvalModPoly.hpp"
#include "basis/poly/operation/common/OperationEvalPoly.hpp"
#include "basis/prewavelet/operation/common/OperationEvalPrewavelet.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalLinear.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationMultipleEvalLinearBoundary.hpp"
#include "basis/modbspline/operation/datadriven/OperationMultipleEvalModBspline.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalModLinear.hpp"
#include "basis/modpoly/operation/datadriven/OperationMultipleEvalModPoly.hpp"
#include "basis/modwavelet/operation/datadriven/OperationMultipleEvalModWavelet.hpp"
#include "basis/poly/operation/datadriven/OperationMultipleEvalPoly.hpp"
#include "basis/prewavelet/operation/datadriven/OperationMultipleEvalPrewavelet.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationTestLinear.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"
#include "basis/modbspline/operation/datadriven/OperationTestModBspline.hpp"
#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"
#include "basis/modpoly/operation/datadriven/OperationTestModPoly.hpp"
#include "basis/modwavelet/operation/datadriven/OperationTestModWavelet.hpp"
#include "basis/poly/operation/datadriven/OperationTestPoly.hpp"
#include "basis/prewavelet/operation/datadriven/OperationTestPrewavelet.hpp"

#include "basis/poly/operation/common/OperationHierarchisationPoly.hpp"
#include "basis/modpoly/operation/common/OperationHierarchisationModPoly.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
#include "basis/linear/noboundary/operation/common/OperationHierarchisationLinear.hpp"
#include "basis/modbspline/operation/common/OperationHierarchisationModBspline.hpp"
#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"
#include "basis/modwavelet/operation/common/OperationHierarchisationModWavelet.hpp"
#include "basis/prewavelet/operation/common/OperationHierarchisationPrewavelet.hpp"

#include "basis/linear/noboundary/operation/common/OperationQuadratureLinear.hpp"
#include "basis/poly/operation/common/OperationQuadraturePoly.hpp"

#include "basis/prewavelet/operation/common/OperationConvertPrewavelet.hpp"

#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLogLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"
#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLogLinear.hpp"
#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLinear.hpp"

#include "basis/linearstretched/boundary/operation/datadriven/OperationTestLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/datadriven/OperationMultipleEvalLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/common/OperationEvalLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/common/OperationHierarchisationLinearStretchedBoundary.hpp"
// @todo (heinecke) removed this when done
#include "basis/linearstretched/boundary/operation/common/OperationUpDownTestLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/operation/pde/OperationLaplaceLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLBLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLDLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLELinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLFLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/OperationLTwoDotProductLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLogLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLogLinearStretchedBoundary.hpp"

// Include all operations on the linearstretched grid
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationTestLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationEvalLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationHierarchisationLinearStretched.hpp"

#include "basis/linearstretched/noboundary/operation/pde/OperationLaplaceLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/OperationLTwoDotProductLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLELinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLBLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLFLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLDLinearStretched.hpp"

#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLogLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLogLinearStretched.hpp"

using namespace sg::base;

namespace sg
{

  namespace GridOperationFactory
  {
#ifdef SG_PDE
    using namespace sg::pde;
    /**
     * gets a pointer to OperationLaplace (OperationMatrix) object
     *
     * @return point to the OperationLaplace object
     */
    //  extern OperationMatrix* createOperationLaplace(Grid& grid);

    /**
     * gets a pointer to OperationLaplace (OperationMatrix) object
     *
     * @return point to the OperationLaplace object
     */
    //  extern OperationMatrix* createOperationLaplace(Grid& grid, sg::base::DataVector& coef);

    /**
     * gets a pointer to OperationLaplace (OperationMatrix) object
     *
     * @return point to the OperationLaplace object
     */
    static OperationMatrix* createOperationLaplace(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLaplaceLinear(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationLaplaceLinearBoundary(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "modlinear") == 0 )
        {
          return new OperationLaplaceModLinear(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "prewavelet") == 0 )
        {
          return new OperationLaplacePrewavelet(grid.getStorage(),
                                                ((PrewaveletGrid*) &grid)->getShadowStorage());
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0 )
        {
          return new OperationLaplaceLinearStretched(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
        {
          return new OperationLaplaceLinearStretchedBoundary(grid.getStorage());
        }
      else
        {
          throw factory_exception("OperationLaplace is not implemented for this grid type.");
        }
    }

    /**
     * gets a pointer to OperationLaplace (OperationMatrix) object
     *
     * @return point to the OperationLaplace object
     */
    static OperationMatrix* createOperationLaplace(Grid& grid, sg::base::DataVector& coef)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLaplaceLinear(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationLaplaceLinearBoundary(grid.getStorage(), coef);
        }
      else
        {
          throw factory_exception("OperationLaplace (with coefficients) is not implemented for this grid type.");
        }
    }

    /**
     * gets a pointer to OperationLTwoDotProduct (OperationMatrix) object
     *
     * @return pointer to OperationLTwoDotProduct object
     */
    static OperationMatrix* createOperationLTwoDotProduct(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLTwoDotProductLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationLTwoDotProductLinearBoundary(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0)
        {
          return new OperationLTwoDotProductLinearStretched(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
        {
          return new OperationLTwoDotProductLinearStretchedBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    static OperationMatrix* createOperationUpDownTest(Grid& grid)
    {
      if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
        {
          return new OperationUpDownTestLinearStretchedBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

#endif //SG_PDE

#ifdef SG_FINANCE
    using namespace sg::finance;
    /**
     * this operation allows you to calculate the following bilinear form
     * needed to solve the multidimensional Black Scholes Equation
     *
     * \f$ \int_{\Omega} S_i S_j \frac{\partial u(\vec{s}}{\partial S_i} \frac{\partial v(\vec{s}}{\partial S_j} d \vec{s}\f$
     *
     * @param coef reference to a DataMatrix object that contains the constant coeffecients of this bilinear from
     */
    static OperationMatrix* createOperationGamma(Grid& grid, DataMatrix& coef)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationGammaLinear(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0)
        {
          return new OperationGammaLinearStretched(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
        {
          return new OperationGammaLinearStretchedBoundary(grid.getStorage(), coef);
        }

      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationGammaLinearBoundary(grid.getStorage(), coef);
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");

    }

    /**
     * this operation allows you to calculate the following bilinear form
     * needed to solve the multidimensional log-transformed Black Scholes Equation
     *
     * \f$ \int_{\Omega} \frac{\partial u(\vec{s}}{\partial S_i} \frac{\partial v(\vec{s}}{\partial S_j} d \vec{s}\f$
     *
     * @param coef reference to a DataVector object that contains the constant coeffecients of this bilinear from
     */
    static OperationMatrix* createOperationGammaLog(Grid& grid, DataMatrix& coef)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationGammaLogLinear(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0)
        {
          return new OperationGammaLogLinearStretched(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
        {
          return new OperationGammaLogLinearStretchedBoundary(grid.getStorage(), coef);
        }

      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationGammaLogLinearBoundary(grid.getStorage(), coef);
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");

    }


    static OperationMatrix* createOperationLB(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLBLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new  OperationLBLinearBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    static OperationMatrix* createOperationLE(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLELinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new  OperationLELinearBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    static OperationMatrix* createOperationLD(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLDLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new  OperationLDLinearBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    static OperationMatrix* createOperationLF(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationLFLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new  OperationLFLinearBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }
    /**
     * this operation allows you to calculate the following bilinear form
     * needed to solve the multidimensional Black Scholes Equation
     *
     * \f$ \int_{\Omega} S_i v(\vec{s}) \frac{\partial u(\vec{s}}{\partial S_i} d \vec{s}\f$
     *
     * @param coef reference to a DataVector object that contains the constant coeffecients of this bilinear from
     */
    static OperationMatrix* createOperationDelta(Grid& grid, DataVector& coef)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationDeltaLinear(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0)
        {
          return new OperationDeltaLinearStretched(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
        {
          return new OperationDeltaLinearStretchedBoundary(grid.getStorage(), coef);
        }

      else if(strcmp(grid.getType(), "linearBoundary") == 0
              || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationDeltaLinearBoundary(grid.getStorage(), coef);
        }

      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    /**
     * this operation allows you to calculate the following bilinear form
     * needed to solve the multidimensional log-transformed Black Scholes Equation
     *
     * \f$ \int_{\Omega} \frac{\partial u(\vec{s}}{\partial S_i} v(\vec{s}) d \vec{s}\f$
     *
     * @param coef reference to a DataVector object that contains the constant coeffecients of this bilinear from
     */
    static OperationMatrix* createOperationDeltaLog(Grid& grid, DataVector& coef)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationDeltaLogLinear(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0)
        {
          return new OperationDeltaLogLinearStretched(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0)
        {
          return new OperationDeltaLogLinearStretchedBoundary(grid.getStorage(), coef);
        }
      else if(strcmp(grid.getType(), "linearBoundary") == 0
              || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0 )
        {
          return new OperationDeltaLogLinearBoundary(grid.getStorage(), coef);
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

#endif //SG_FINANCE



    /**
     * gets a pointer to OperationEval object
     *
     * @return pointer to the OperationEval object
     */
    static OperationEval* createOperationEval(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationEvalLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0
              || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
              || strcmp(grid.getType(), "TruncatedTrapezoid") == 0
              || strcmp(grid.getType(), "squareRoot") == 0)
        {
          return new OperationEvalLinearBoundary(grid.getStorage());
        }

      else if(strcmp(grid.getType(), "modBspline") == 0 )
        {
          return new OperationEvalModBspline(grid.getStorage(),
                                             ((ModBsplineGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modlinear") == 0 )
        {
          return new OperationEvalModLinear(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "modpoly") == 0 )
        {
          return new OperationEvalModPoly(grid.getStorage(),
                                          ((ModPolyGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modWavelet") == 0 )
        {
          return new OperationEvalModWavelet(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "poly") == 0 )
        {
          return new OperationEvalPoly(grid.getStorage(),
                                       ((PolyGrid*) &grid)->getDegree());
        }

      else if(strcmp(grid.getType(), "prewavelet") == 0 )
        {
          return new OperationEvalPrewavelet(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0 )
        {
          return new OperationEvalLinearStretched(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
        {
          return new OperationEvalLinearStretchedBoundary(grid.getStorage());
        }
      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    /**
     * gets a pointer to OperationMultipleEval object
     *
     * @param dataset the dataset that should be evaluated on the sparse grid
     *
     * @return pointer to the OperationMultipleEval object
     */
    static OperationMultipleEval* createOperationMultipleEval(Grid& grid, DataMatrix* dataset)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationMultipleEvalLinear(grid.getStorage(), dataset);
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0
              || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationMultipleEvalLinearBoundary(grid.getStorage(), dataset);
        }

      else if(strcmp(grid.getType(), "modBspline") == 0 )
        {
          return new OperationMultipleEvalModBspline(grid.getStorage(),
                                                     ((ModBsplineGrid*) &grid)->getDegree(), dataset);
        }
      else if(strcmp(grid.getType(), "modlinear") == 0 )
        {
          return new OperationMultipleEvalModLinear(grid.getStorage(), dataset);
        }
      else if(strcmp(grid.getType(), "modpoly") == 0 )
        {
          return new OperationMultipleEvalModPoly(grid.getStorage(),
                                                  ((ModPolyGrid*) &grid)->getDegree(), dataset);
        }
      else if(strcmp(grid.getType(), "modWavelet") == 0 )
        {
          return new OperationMultipleEvalModWavelet(grid.getStorage(), dataset);
        }
      else if(strcmp(grid.getType(), "poly") == 0 )
        {
          return new OperationMultipleEvalPoly(grid.getStorage(),
                                               ((PolyGrid*) &grid)->getDegree(), dataset);
        }
      else if(strcmp(grid.getType(), "prewavelet") == 0 )
        {
          return new OperationMultipleEvalPrewavelet(grid.getStorage(), dataset);
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0 )
        {
          return new OperationMultipleEvalLinearStretched(grid.getStorage(), dataset);
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
        {
          return new OperationMultipleEvalLinearStretchedBoundary(grid.getStorage(), dataset);
        }

      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }


#ifdef SG_DATADRIVEN
    using namespace sg::datadriven;
    /**
     * gets a pointer to OperationTest object
     *
     * @return pointer to the OperationTest object
     */
    static OperationTest* createOperationTest(Grid& grid)
    {
      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationTestLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0
              || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
        {
          return new OperationTestLinearBoundary(grid.getStorage());
        }

      else if(strcmp(grid.getType(), "modBspline") == 0 )
        {
          return new OperationTestModBspline(grid.getStorage(),
                                             ((ModBsplineGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modlinear") == 0 )
        {
          return new OperationTestModLinear(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "modpoly") == 0 )
        {
          return new OperationTestModPoly(grid.getStorage(),
                                          ((ModPolyGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modWavelet") == 0 )
        {
          return new OperationTestModWavelet(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "poly") == 0 )
        {
          return new OperationTestPoly(grid.getStorage(),
                                       ((PolyGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "prewavelet") == 0 )
        {
          return new OperationTestPrewavelet(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0 )
        {
          return new OperationTestLinearStretched(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
        {
          return new OperationTestLinearStretchedBoundary(grid.getStorage());
        }

      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }
#endif //SG_DATADRIVEN

    /**
     * gets a pointer to OperationHierarchisation object
     *
     * @return pointer to the OperationHierarchisation object
     */
    static OperationHierarchisation* createOperationHierarchisation(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationHierarchisationLinear(grid.getStorage());
        }


      else if(strcmp(grid.getType(), "linearBoundary") == 0
              || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
              || strcmp(grid.getType(), "TruncatedTrapezoid") == 0
              || strcmp(grid.getType(), "squareRoot") == 0)
        {
          return new OperationHierarchisationLinearBoundary(grid.getStorage());
        }

      else if(strcmp(grid.getType(), "modBspline") == 0 )
        {
          return new OperationHierarchisationModBspline(grid.getStorage(),
                                                        ((ModBsplineGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modlinear") == 0 )
        {
          return new OperationHierarchisationModLinear(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "poly") == 0 )
        {
          return new OperationHierarchisationPoly(grid.getStorage(),
                                                  ((PolyGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modpoly") == 0 )
        {
          return new OperationHierarchisationModPoly(grid.getStorage(),
                                                     ((ModPolyGrid*) &grid)->getDegree());
        }
      else if(strcmp(grid.getType(), "modWavelet") == 0 )
        {
          return new OperationHierarchisationModWavelet(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretched") == 0 )
        {
          return new OperationHierarchisationLinearStretched(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "linearStretchedTrapezoidBoundary") == 0 )
        {
          return new OperationHierarchisationLinearStretchedBoundary(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "prewavelet") == 0 )
        {
          return new OperationHierarchisationPrewavelet(grid.getStorage(),
                                                        ((PrewaveletGrid*) &grid)->getShadowStorage());
        }

      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }

    /**
     * Creates a OperationQuadrature
     *
     * @return pointer to the OperationQuadrature object
     */
    static OperationQuadrature* createOperationQuadrature(Grid& grid)
    {

      if(strcmp(grid.getType(), "linear") == 0)
        {
          return new OperationQuadratureLinear(grid.getStorage());
        }
      else if(strcmp(grid.getType(), "poly") == 0 )
        {
          if(((PolyGrid*) &grid)->getDegree()>3) {
            throw factory_exception("OperationQuadrature is not implemented for polynomials with degree higher than 3.");
          }
          else {
            return new OperationQuadraturePoly(grid.getStorage(), ((PolyGrid*) &grid)->getDegree());
          }
        }
      else
        throw factory_exception("OperationQuadrature is not implemented for this grid type.");
    }

    /**
     * gets a pointer to OperationConvert object
     *
     * @return pointer to the OperationConvert object
     */
    static OperationConvert* createOperationConvert(Grid& grid)
    {
      if(strcmp(grid.getType(), "prewavelet") == 0 )
        {
          return new OperationConvertPrewavelet(grid.getStorage(),
                                                ((PrewaveletGrid*) &grid)->getShadowStorage());
        }

      else
        throw factory_exception("OperationLaplace is not implemented for this grid type.");
    }


    /**
     * gets a pointer to OperationIdentity (OperationMatrix) object
     *
     * @return point to the OperationIdentity object
     */
    static OperationMatrix* createOperationIdentity(Grid& grid)
    {
      return new OperationIdentity();
    }

  }

}

#endif /*OPERATIONS_FACTORY_HPP*/
