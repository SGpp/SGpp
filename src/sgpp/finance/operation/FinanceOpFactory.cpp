/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "finance/operation/FinanceOpFactory.hpp"

#include <cstring>

#include "base/exception/factory_exception.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationGammaLinear.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationGammaLinearStretched.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationGammaLogLinearStretchedBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationGammaLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationGammaLogLinear.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationGammaLogLinearStretched.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationGammaLinearStretchedBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationGammaLogLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationLBLinear.hpp"
#include "finance/basis/linear/boundary/operation/OperationLBLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationLELinear.hpp"
#include "finance/basis/linear/boundary/operation/OperationLELinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationLDLinear.hpp"
#include "finance/basis/linear/boundary/operation/OperationLDLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationLFLinear.hpp"
#include "finance/basis/linear/boundary/operation/OperationLFLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationDeltaLinear.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationDeltaLinearStretched.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationDeltaLinearStretchedBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationDeltaLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationDeltaLogLinear.hpp"
#include "finance/basis/linearstretched/noboundary/operation/OperationDeltaLogLinearStretched.hpp"
#include "finance/basis/linearstretched/boundary/operation/OperationDeltaLogLinearStretchedBoundary.hpp"
#include "finance/basis/linear/boundary/operation/OperationDeltaLogLinearBoundary.hpp"

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonBLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonCLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonDLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonELinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonFLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonGLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonHLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonKLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonWLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonXLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonYLinear.hpp"
#include "finance/basis/linear/noboundary/operation/pde/OperationHestonZLinear.hpp"

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

base::OperationMatrix* createOperationHestonBLog(base::Grid& grid, base::DataMatrix& coef)
{
	return new finance::OperationHestonBLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonCLog(base::Grid& grid, base::DataMatrix& coef)
{
	return new finance::OperationHestonCLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonDLog(base::Grid& grid, base::DataVector& coef)
{
	return new finance::OperationHestonDLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonELog(base::Grid& grid, base::DataVector& coef)
{
	return new finance::OperationHestonELinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonFLog(base::Grid& grid, base::DataVector& coef)
{
	return new finance::OperationHestonFLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonGLog(base::Grid& grid, base::DataVector& coef)
{
	return new finance::OperationHestonGLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonHLog(base::Grid& grid, base::DataMatrix& coef)
{
	return new finance::OperationHestonHLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonKLog(base::Grid& grid, double***** coef)
{
	return new finance::OperationHestonKLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonX(base::Grid& grid, base::DataMatrix& coef)
{
	return new finance::OperationHestonXLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonY(base::Grid& grid, base::DataMatrix& coef)
{
	return new finance::OperationHestonYLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonW(base::Grid& grid, base::DataMatrix& coef)
{
	return new finance::OperationHestonWLinear(grid.getStorage(), coef);
}

base::OperationMatrix* createOperationHestonZ(base::Grid& grid, base::DataVector& coef)
{
	return new finance::OperationHestonZLinear(grid.getStorage(), coef);
}

}
}

