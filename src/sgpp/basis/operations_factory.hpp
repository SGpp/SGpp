/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

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


#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSSELinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPSSELinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeAVXLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPAVXLinear.hpp"

#ifdef USEOCL
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPHybridSSEOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeHybridSSEOCLLinear.hpp"
#endif

#ifdef USEARBB
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeArBBLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#endif


#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
#include "basis/linear/noboundary/operation/common/OperationHierarchisationLinear.hpp"
#include "basis/modbspline/operation/common/OperationHierarchisationModBspline.hpp"
#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"
#include "basis/modpoly/operation/common/OperationHierarchisationModPoly.hpp"
#include "basis/modwavelet/operation/common/OperationHierarchisationModWavelet.hpp"
#include "basis/poly/operation/common/OperationHierarchisationPoly.hpp"
#include "basis/prewavelet/operation/common/OperationHierarchisationPrewavelet.hpp"

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
	static OperationMatrix* createOperationLaplace(Grid& grid_type)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationLaplaceLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationLaplaceLinearBoundary(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "modlinear") == 0 )
		{
			return new OperationLaplaceModLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "prewavelet") == 0 )
		{
			return new OperationLaplacePrewavelet(grid_type.getStorage(),
					((PrewaveletGrid*) &grid_type)->getShadowStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0 )
		{
			return new OperationLaplaceLinearStretched(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0 )
		{
			return new OperationLaplaceLinearStretchedBoundary(grid_type.getStorage());
		}

		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");


	}

	/**
	 * gets a pointer to OperationLTwoDotProduct (OperationMatrix) object
	 *
	 * @return pointer to OperationLTwoDotProduct object
	 */
	static OperationMatrix* createOperationLTwoDotProduct(Grid& grid_type)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationLTwoDotProductLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationLTwoDotProductLinearBoundary(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0)
		{
			return new OperationLTwoDotProductLinearStretched(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0)
		{
			return new OperationLTwoDotProductLinearStretchedBoundary(grid_type.getStorage());
		}
		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}

	static OperationMatrix* createOperationUpDownTest(Grid& grid_type)
	{
		if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0)
		{
			return new OperationUpDownTestLinearStretchedBoundary(grid_type.getStorage());
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
	static OperationMatrix* createOperationGamma(Grid& grid_type, DataMatrix& coef)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationGammaLinear(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0)
		{
			return new OperationGammaLinearStretched(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0)
		{
			return new OperationGammaLinearStretchedBoundary(grid_type.getStorage(), coef);
		}

		else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationGammaLinearBoundary(grid_type.getStorage(), coef);
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
	static OperationMatrix* createOperationGammaLog(Grid& grid_type, DataMatrix& coef)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationGammaLogLinear(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0)
		{
			return new OperationGammaLogLinearStretched(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0)
		{
			return new OperationGammaLogLinearStretchedBoundary(grid_type.getStorage(), coef);
		}

		else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationGammaLogLinearBoundary(grid_type.getStorage(), coef);
		}
		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");

	}


	static OperationMatrix* createOperationLB(Grid& grid_type)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationLBLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new  OperationLBLinearBoundary(grid_type.getStorage());
		}
		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}

	static OperationMatrix* createOperationLE(Grid& grid_type)
		{

			if(strcmp(grid_type.getType(), "linear") == 0)
			{
				return new OperationLELinear(grid_type.getStorage());
			}


			else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
			{
				return new  OperationLELinearBoundary(grid_type.getStorage());
			}
			else
				throw factory_exception("OperationLaplace is not implemented for this grid type.");
		}

	static OperationMatrix* createOperationLD(Grid& grid_type)
		{

			if(strcmp(grid_type.getType(), "linear") == 0)
			{
				return new OperationLDLinear(grid_type.getStorage());
			}


			else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
			{
				return new  OperationLDLinearBoundary(grid_type.getStorage());
			}
			else
				throw factory_exception("OperationLaplace is not implemented for this grid type.");
		}

	static OperationMatrix* createOperationLF(Grid& grid_type)
		{

			if(strcmp(grid_type.getType(), "linear") == 0)
			{
				return new OperationLFLinear(grid_type.getStorage());
			}


			else if(strcmp(grid_type.getType(), "linearBoundary") == 0 || strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
			{
				return new  OperationLFLinearBoundary(grid_type.getStorage());
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
	static OperationMatrix* createOperationDelta(Grid& grid_type, DataVector& coef)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationDeltaLinear(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0)
		{
			return new OperationDeltaLinearStretched(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0)
		{
			return new OperationDeltaLinearStretchedBoundary(grid_type.getStorage(), coef);
		}

		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationDeltaLinearBoundary(grid_type.getStorage(), coef);
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
	static OperationMatrix* createOperationDeltaLog(Grid& grid_type, DataVector& coef)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationDeltaLogLinear(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0)
		{
			return new OperationDeltaLogLinearStretched(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0)
		{
			return new OperationDeltaLogLinearStretchedBoundary(grid_type.getStorage(), coef);
		}
		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0 )
		{
			return new OperationDeltaLogLinearBoundary(grid_type.getStorage(), coef);
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
	static OperationEval* createOperationEval(Grid& grid_type)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationEvalLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0
				|| strcmp(grid_type.getType(), "TruncatedTrapezoid") == 0
				|| strcmp(grid_type.getType(), "squareRoot") == 0)
		{
			return new OperationEvalLinearBoundary(grid_type.getStorage());
		}

		else if(strcmp(grid_type.getType(), "modBspline") == 0 )
		{
			return new OperationEvalModBspline(grid_type.getStorage(),
					((ModBsplineGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "modlinear") == 0 )
		{
			return new OperationEvalModLinear(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "modpoly") == 0 )
		{
			return new OperationEvalModPoly(grid_type.getStorage(),
					((ModPolyGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "modWavelet") == 0 )
		{
			return new OperationEvalModWavelet(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "poly") == 0 )
		{
			return new OperationEvalPoly(grid_type.getStorage(),
					((PolyGrid*) &grid_type)->getDegree());
		}

		else if(strcmp(grid_type.getType(), "prewavelet") == 0 )
		{
			return new OperationEvalPrewavelet(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0 )
		{
			return new OperationEvalLinearStretched(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0 )
		{
			return new OperationEvalLinearStretchedBoundary(grid_type.getStorage());
		}
		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}

	/**
	 * gets a pointer to OperationMultipleEval object
	 *
	 * @param dataset the dataset that should be evaluated on the sparse grid
	 *
	 * @return pointer to the OperationB object
	 */
	static OperationMultipleEval* createOperationMultipleEval(Grid& grid_type, DataMatrix* dataset)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationMultipleEvalLinear(grid_type.getStorage(), dataset);
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationMultipleEvalLinearBoundary(grid_type.getStorage(), dataset);
		}

		else if(strcmp(grid_type.getType(), "modBspline") == 0 )
		{
			return new OperationMultipleEvalModBspline(grid_type.getStorage(),
					((ModBsplineGrid*) &grid_type)->getDegree(), dataset);
		}
		else if(strcmp(grid_type.getType(), "modlinear") == 0 )
		{
			return new OperationMultipleEvalModLinear(grid_type.getStorage(), dataset);
		}
		else if(strcmp(grid_type.getType(), "modpoly") == 0 )
		{
			return new OperationMultipleEvalModPoly(grid_type.getStorage(),
					((ModPolyGrid*) &grid_type)->getDegree(), dataset);
		}
		else if(strcmp(grid_type.getType(), "modWavelet") == 0 )
		{
			return new OperationMultipleEvalModWavelet(grid_type.getStorage(), dataset);
		}
		else if(strcmp(grid_type.getType(), "poly") == 0 )
		{
			return new OperationMultipleEvalPoly(grid_type.getStorage(),
					((PolyGrid*) &grid_type)->getDegree(), dataset);
		}
		else if(strcmp(grid_type.getType(), "prewavelet") == 0 )
		{
			return new OperationMultipleEvalPrewavelet(grid_type.getStorage(), dataset);
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0 )
		{
			return new OperationMultipleEvalLinearStretched(grid_type.getStorage(), dataset);
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0 )
		{
			return new OperationMultipleEvalLinearStretchedBoundary(grid_type.getStorage(), dataset);
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
	static OperationTest* createOperationTest(Grid& grid_type)
	{
		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationTestLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			return new OperationTestLinearBoundary(grid_type.getStorage());
		}

		else if(strcmp(grid_type.getType(), "modBspline") == 0 )
		{
			return new OperationTestModBspline(grid_type.getStorage(),
					((ModBsplineGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "modlinear") == 0 )
		{
			return new OperationTestModLinear(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "modpoly") == 0 )
		{
			return new OperationTestModPoly(grid_type.getStorage(),
					((ModPolyGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "modWavelet") == 0 )
		{
			return new OperationTestModWavelet(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "poly") == 0 )
		{
			return new OperationTestPoly(grid_type.getStorage(),
					((PolyGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "prewavelet") == 0 )
		{
			return new OperationTestPrewavelet(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0 )
		{
			return new OperationTestLinearStretched(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0 )
		{
			return new OperationTestLinearStretchedBoundary(grid_type.getStorage());
		}

		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}
#endif //SG_DATADRIVEN



#ifdef SG_PARALLEL
//using namespace sg::parallel;
	/**
	 * gets a pointer to OperationBVectorized object
	 *
	 * @param VecType Type of Vectorization used: Currently supported: SSE, AVX
	 * @param dataset the dataset that should be evaluated on the sparse grid
	 *
	 * @return pointer to the OperationB object
	 */
	static OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(Grid& grid_type, const std::string& VecType, DataMatrix* dataset)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			if (VecType == "SSE")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSSELinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "AVX")
			{
				return new sg::parallel::OperationMultipleEvalIterativeAVXLinear(grid_type.getStorage(), dataset);
			}
		#ifdef USEOCL
			else if (VecType == "OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeOCLLinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "HYBRID_SSE_OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeHybridSSEOCLLinear(grid_type.getStorage(), dataset);
			}
		#endif
		#ifdef USEARBB
			else if (VecType == "ArBB")
			{
				return new sg::parallel::OperationMultipleEvalIterativeArBBLinear(grid_type.getStorage(), dataset);
			}
		#endif
			else
			{
				throw factory_exception("Unsupported vectorization type");
			}
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			if (VecType == "SSE")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSSELinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "AVX")
			{
				return new sg::parallel::OperationMultipleEvalIterativeAVXLinear(grid_type.getStorage(), dataset);
			}
		#ifdef USEOCL
			else if (VecType == "OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeOCLLinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "HYBRID_SSE_OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeHybridSSEOCLLinear(grid_type.getStorage(), dataset);
			}
		#endif
		#ifdef USEARBB
			else if (VecType == "ArBB")
			{
				return new sg::parallel::OperationMultipleEvalIterativeArBBLinear(grid_type.getStorage(), dataset);
			}
		#endif
			else
			{
				throw factory_exception("Unsupported vectorization type");
			}
		}

		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}

	/**
	 * gets a pointer to OperationBVectorizedSP object
	 *
	 * @param VecType Type of Vectorization used: Currently supported: SSE, AVX
	 * @param dataset the dataset that should be evaluated on the sparse grid
	 *
	 * @return pointer to the OperationBSP object
	 */
	static OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(Grid& grid_type, const std::string& VecType, DataMatrixSP* dataset)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			if (VecType == "SSE")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPSSELinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "AVX")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPAVXLinear(grid_type.getStorage(), dataset);
			}
		#ifdef USEOCL
			else if (VecType == "OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPOCLLinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "HYBRID_SSE_OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPHybridSSEOCLLinear(grid_type.getStorage(), dataset);
			}
		#endif
		#ifdef USEARBB
			else if (VecType == "ArBB")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPArBBLinear(grid_type.getStorage(), dataset);
			}
		#endif
			else
			{
				throw factory_exception("Unsupported vectorization type");
			}
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0)
		{
			if (VecType == "SSE")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPSSELinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "AVX")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPAVXLinear(grid_type.getStorage(), dataset);
			}
		#ifdef USEOCL
			else if (VecType == "OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPOCLLinear(grid_type.getStorage(), dataset);
			}
			else if (VecType == "HYBRID_SSE_OCL")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPHybridSSEOCLLinear(grid_type.getStorage(), dataset);
			}
		#endif
		#ifdef USEARBB
			else if (VecType == "ArBB")
			{
				return new sg::parallel::OperationMultipleEvalIterativeSPArBBLinear(grid_type.getStorage(), dataset);
			}
		#endif
			else
			{
				throw factory_exception("Unsupported vectorization type");
			}
		}

		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}
#endif //SG_PARALLEL

	/**
	 * gets a pointer to OperationHierarchisation object
	 *
	 * @return pointer to the OperationHierarchisation object
	 */
	static OperationHierarchisation* createOperationHierarchisation(Grid& grid_type)
	{

		if(strcmp(grid_type.getType(), "linear") == 0)
		{
			return new OperationHierarchisationLinear(grid_type.getStorage());
		}


		else if(strcmp(grid_type.getType(), "linearBoundary") == 0
				|| strcmp(grid_type.getType(), "linearTrapezoidBoundary") == 0
				|| strcmp(grid_type.getType(), "TruncatedTrapezoid") == 0
				|| strcmp(grid_type.getType(), "squareRoot") == 0)
		{
			return new OperationHierarchisationLinearBoundary(grid_type.getStorage());
		}

		else if(strcmp(grid_type.getType(), "modBspline") == 0 )
		{
			return new OperationHierarchisationModBspline(grid_type.getStorage(),
					((ModBsplineGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "modlinear") == 0 )
		{
			return new OperationHierarchisationModLinear(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "modpoly") == 0 )
		{
			return new OperationHierarchisationModPoly(grid_type.getStorage(),
					((ModPolyGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "modWavelet") == 0 )
		{
			return new OperationHierarchisationModWavelet(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretched") == 0 )
		{
			return new OperationHierarchisationLinearStretched(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "linearStretchedTrapezoidBoundary") == 0 )
		{
			return new OperationHierarchisationLinearStretchedBoundary(grid_type.getStorage());
		}
		else if(strcmp(grid_type.getType(), "poly") == 0 )
		{
			return new OperationHierarchisationPoly(grid_type.getStorage(),
					((PolyGrid*) &grid_type)->getDegree());
		}
		else if(strcmp(grid_type.getType(), "prewavelet") == 0 )
		{
			return new OperationHierarchisationPrewavelet(grid_type.getStorage(),
					((PrewaveletGrid*) &grid_type)->getShadowStorage());
		}

		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}

	/**
	 * gets a pointer to OperationConvert object
	 *
	 * @return pointer to the OperationConvert object
	 */
	static OperationConvert* createOperationConvert(Grid& grid_type)
	{
		if(strcmp(grid_type.getType(), "prewavelet") == 0 )
		{
			return new OperationConvertPrewavelet(grid_type.getStorage(),
					((PrewaveletGrid*) &grid_type)->getShadowStorage());
		}

		else
			throw factory_exception("OperationLaplace is not implemented for this grid type.");
	}


	/**
	 * gets a pointer to OperationIdentity (OperationMatrix) object
	 *
	 * @return point to the OperationIdentity object
	 */
	static OperationMatrix* createOperationIdentity(Grid& grid_type)
	{
		return new OperationIdentity();
	}

}

}

#endif /*OPERATIONS_FACTORY_HPP*/
