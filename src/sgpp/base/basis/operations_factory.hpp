/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONS_FACTORY_HPP
#define OPERATIONS_FACTORY_HPP

#include <cstring>

#include "base/exception/factory_exception.hpp"

#include "base/grid/Grid.hpp"
#include "base/grid/type/PrewaveletGrid.hpp"
#include "base/grid/type/PolyGrid.hpp"
#include "base/grid/type/ModPolyGrid.hpp"
#include "base/grid/type/ModBsplineGrid.hpp"

#include "base/operation/OperationMatrix.hpp"
#include "base/operation/OperationEval.hpp"
#include "datadriven/operation/OperationMultipleEval.hpp"
#include "pde/operation/OperationParabolicPDESolverSystem.hpp"
#include "datadriven/operation/OperationTest.hpp"
#include "datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "base/operation/OperationHierarchisation.hpp"
#include "base/operation/OperationQuadrature.hpp"
#include "base/operation/OperationConvert.hpp"
#include "base/operation/OperationIdentity.hpp"

#include "pde/basis/linear/noboundary/operation/OperationLaplaceLinear.hpp"
#include "pde/basis/linear/boundary/operation/OperationLaplaceLinearBoundary.hpp"
#include "pde/basis/modlinear/operation/OperationLaplaceModLinear.hpp"
#include "datadriven/basis/prewavelet/operation/OperationLaplacePrewavelet.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationGammaLinear.hpp"
#include "finance/basis/linear/boundary/operation/OperationGammaLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationGammaLogLinear.hpp"
#include "finance/basis/linear/boundary/operation/OperationGammaLogLinearBoundary.hpp"

#include "finance/basis/linear/boundary/operation/OperationLELinearBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationLBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationLFLinearBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationLDLinearBoundary.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationLELinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationLBLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationLFLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationLDLinear.hpp"

#include "pde/basis/linear/noboundary/operation/OperationLTwoDotProductLinear.hpp"
#include "pde/basis/linear/boundary/operation/OperationLTwoDotProductLinearBoundary.hpp"

#include "base/basis/modwavelet/operation/OperationEvalModWavelet.hpp"
#include "base/basis/linear/noboundary/operation/OperationEvalLinear.hpp"
#include "base/basis/linear/boundary/operation/OperationEvalLinearBoundary.hpp"
#include "base/basis/modbspline/operation/OperationEvalModBspline.hpp"
#include "base/basis/modlinear/operation/OperationEvalModLinear.hpp"
#include "base/basis/modpoly/operation/OperationEvalModPoly.hpp"
#include "base/basis/poly/operation/OperationEvalPoly.hpp"
#include "base/basis/prewavelet/operation/OperationEvalPrewavelet.hpp"

#include "datadriven/basis/linear/noboundary/operation/OperationMultipleEvalLinear.hpp"
#include "datadriven/basis/linear/boundary/operation/OperationMultipleEvalLinearBoundary.hpp"
#include "datadriven/basis/modbspline/operation/OperationMultipleEvalModBspline.hpp"
#include "datadriven/basis/modlinear/operation/OperationMultipleEvalModLinear.hpp"
#include "datadriven/basis/modpoly/operation/OperationMultipleEvalModPoly.hpp"
#include "datadriven/basis/modwavelet/operation/OperationMultipleEvalModWavelet.hpp"
#include "datadriven/basis/poly/operation/OperationMultipleEvalPoly.hpp"
#include "datadriven/basis/prewavelet/operation/OperationMultipleEvalPrewavelet.hpp"

#include "datadriven/basis/linear/noboundary/operation/OperationTestLinear.hpp"
#include "datadriven/basis/linear/boundary/operation/OperationTestLinearBoundary.hpp"
#include "datadriven/basis/modbspline/operation/OperationTestModBspline.hpp"
#include "datadriven/basis/modlinear/operation/OperationTestModLinear.hpp"
#include "datadriven/basis/modpoly/operation/OperationTestModPoly.hpp"
#include "datadriven/basis/modwavelet/operation/OperationTestModWavelet.hpp"
#include "datadriven/basis/poly/operation/OperationTestPoly.hpp"
#include "datadriven/basis/prewavelet/operation/OperationTestPrewavelet.hpp"

#include "base/basis/poly/operation/OperationHierarchisationPoly.hpp"
#include "base/basis/modpoly/operation/OperationHierarchisationModPoly.hpp"
#include "base/basis/linear/boundary/operation/OperationHierarchisationLinearBoundary.hpp"
#include "base/basis/linear/noboundary/operation/OperationHierarchisationLinear.hpp"
#include "base/basis/modbspline/operation/OperationHierarchisationModBspline.hpp"
#include "base/basis/modlinear/operation/OperationHierarchisationModLinear.hpp"
#include "base/basis/modwavelet/operation/OperationHierarchisationModWavelet.hpp"
#include "base/basis/prewavelet/operation/OperationHierarchisationPrewavelet.hpp"

#include "base/basis/linear/noboundary/operation/OperationQuadratureLinear.hpp"
#include "base/basis/poly/operation/OperationQuadraturePoly.hpp"

#include "base/basis/prewavelet/operation/OperationConvertPrewavelet.hpp"

#include "finance/basis/linear/boundary/operation/OperationDeltaLogLinearBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationDeltaLinearBoundary.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationDeltaLogLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationDeltaLinear.hpp"

#include "datadriven/basis/linearstretched/boundary/operation/OperationTestLinearStretchedBoundary.hpp"
#include "datadriven/basis/linearstretched/boundary/operation/OperationMultipleEvalLinearStretchedBoundary.hpp"
#include "base/basis/linearstretched/boundary/operation/OperationEvalLinearStretchedBoundary.hpp"
#include "base/basis/linearstretched/boundary/operation/OperationHierarchisationLinearStretchedBoundary.hpp"
// @todo (heinecke) removed this when done
#include "base/basis/linearstretched/boundary/operation/OperationUpDownTestLinearStretchedBoundary.hpp"

#include "pde/basis/linearstretched/boundary/operation/OperationLaplaceLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLBLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLDLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLELinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLFLinearStretchedBoundary.hpp"
#include "pde/basis/linearstretched/boundary/operation/OperationLTwoDotProductLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationDeltaLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationGammaLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationDeltaLogLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationGammaLogLinearStretchedBoundary.hpp"

// Include all operations on the linearstretched grid
#include "datadriven/basis/linearstretched/noboundary/operation/OperationMultipleEvalLinearStretched.hpp"
#include "datadriven/basis/linearstretched/noboundary/operation/OperationTestLinearStretched.hpp"
#include "base/basis/linearstretched/noboundary/operation/OperationEvalLinearStretched.hpp"
#include "base/basis/linearstretched/noboundary/operation/OperationHierarchisationLinearStretched.hpp"

#include "pde/basis/linearstretched/noboundary/operation/OperationLaplaceLinearStretched.hpp"
#include "pde/basis/linearstretched/noboundary/operation/OperationLTwoDotProductLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLELinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLBLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLFLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLDLinearStretched.hpp"

#include "finance/basis/linearstretched/noboundary/operation/OperationDeltaLinearStretched.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationGammaLinearStretched.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationDeltaLogLinearStretched.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationGammaLogLinearStretched.hpp"

#include "datadriven/basis/linear/boundary/operation/OperationRegularizationDiagonalLinearBoundary.hpp"


using namespace sg::base;

namespace sg
{

  namespace op_factory //GridOperationFactory
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

  static OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid, int mode, double k)
  {
    if(strcmp(grid.getType(), "linear") == 0
       || strcmp(grid.getType(), "linearBoundary") == 0 
       || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0
       || strcmp(grid.getType(), "modlinear") == 0) {
      return new datadriven::OperationRegularizationDiagonalLinearBoundary(grid.getStorage(), mode, k);
    }
    else
      throw base::factory_exception("OperationRegularizationDiagonal is not implemented for this grid type.");
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
        throw factory_exception("OperationConvert is not implemented for this grid type.");
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
