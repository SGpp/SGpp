/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include <cstring>

#include "opt/operation/OpFactory.hpp"
#include "base/exception/factory_exception.hpp"

#include "opt/grid/LinearGrid.hpp"
#include "opt/grid/LinearTrapezoidBoundaryGrid.hpp"
#include "opt/grid/LinearClenshawCurtisGrid.hpp"
#include "opt/grid/ModLinearGrid.hpp"
#include "opt/grid/BsplineGrid.hpp"
#include "opt/grid/BsplineTrapezoidBoundaryGrid.hpp"
#include "opt/grid/BsplineClenshawCurtisGrid.hpp"
#include "opt/grid/ModBsplineGrid.hpp"
#include "opt/grid/WaveletGrid.hpp"
#include "opt/grid/WaveletTrapezoidBoundaryGrid.hpp"
#include "opt/grid/ModWaveletGrid.hpp"

#include "opt/basis/bspline/boundary/operation/OperationEvalBsplineBoundary.hpp"
#include "opt/basis/bspline/clenshawcurtis/operation/OperationEvalBsplineClenshawCurtis.hpp"
#include "opt/basis/bspline/modified/operation/OperationEvalModBspline.hpp"
#include "opt/basis/bspline/noboundary/operation/OperationEvalBspline.hpp"
#include "opt/basis/linear/boundary/operation/OperationEvalLinearBoundary.hpp"
#include "opt/basis/linear/clenshawcurtis/operation/OperationEvalLinearClenshawCurtis.hpp"
#include "opt/basis/linear/modified/operation/OperationEvalModLinear.hpp"
#include "opt/basis/linear/noboundary/operation/OperationEvalLinear.hpp"
#include "opt/basis/wavelet/boundary/operation/OperationEvalWaveletBoundary.hpp"
#include "opt/basis/wavelet/modified/operation/OperationEvalModWavelet.hpp"
#include "opt/basis/wavelet/noboundary/operation/OperationEvalWavelet.hpp"

#include "opt/basis/bspline/boundary/operation/OperationEvalGradientBsplineBoundary.hpp"
#include "opt/basis/bspline/clenshawcurtis/operation/OperationEvalGradientBsplineClenshawCurtis.hpp"
#include "opt/basis/bspline/modified/operation/OperationEvalGradientModBspline.hpp"
#include "opt/basis/bspline/noboundary/operation/OperationEvalGradientBspline.hpp"
#include "opt/basis/wavelet/boundary/operation/OperationEvalGradientWaveletBoundary.hpp"
#include "opt/basis/wavelet/modified/operation/OperationEvalGradientModWavelet.hpp"
#include "opt/basis/wavelet/noboundary/operation/OperationEvalGradientWavelet.hpp"

#include "opt/basis/bspline/boundary/operation/OperationEvalHessianBsplineBoundary.hpp"
#include "opt/basis/bspline/clenshawcurtis/operation/OperationEvalHessianBsplineClenshawCurtis.hpp"
#include "opt/basis/bspline/modified/operation/OperationEvalHessianModBspline.hpp"
#include "opt/basis/bspline/noboundary/operation/OperationEvalHessianBspline.hpp"
#include "opt/basis/wavelet/boundary/operation/OperationEvalHessianWaveletBoundary.hpp"
#include "opt/basis/wavelet/modified/operation/OperationEvalHessianModWavelet.hpp"
#include "opt/basis/wavelet/noboundary/operation/OperationEvalHessianWavelet.hpp"

#include "opt/basis/bspline/boundary/operation/OperationMultipleHierarchisationBsplineBoundary.hpp"
#include "opt/basis/bspline/clenshawcurtis/operation/OperationMultipleHierarchisationBsplineClenshawCurtis.hpp"
#include "opt/basis/bspline/modified/operation/OperationMultipleHierarchisationModBspline.hpp"
#include "opt/basis/bspline/noboundary/operation/OperationMultipleHierarchisationBspline.hpp"
#include "opt/basis/linear/boundary/operation/OperationMultipleHierarchisationLinearBoundary.hpp"
#include "opt/basis/linear/clenshawcurtis/operation/OperationMultipleHierarchisationLinearClenshawCurtis.hpp"
#include "opt/basis/linear/modified/operation/OperationMultipleHierarchisationModLinear.hpp"
#include "opt/basis/linear/noboundary/operation/OperationMultipleHierarchisationLinear.hpp"
#include "opt/basis/wavelet/boundary/operation/OperationMultipleHierarchisationWaveletBoundary.hpp"
#include "opt/basis/wavelet/modified/operation/OperationMultipleHierarchisationModWavelet.hpp"
#include "opt/basis/wavelet/noboundary/operation/OperationMultipleHierarchisationWavelet.hpp"

#include "opt/basis/bspline/boundary/operation/OperationEvalPartialDerivativeBsplineBoundary.hpp"
#include "opt/basis/bspline/clenshawcurtis/operation/OperationEvalPartialDerivativeBsplineClenshawCurtis.hpp"
#include "opt/basis/bspline/modified/operation/OperationEvalPartialDerivativeModBspline.hpp"
#include "opt/basis/bspline/noboundary/operation/OperationEvalPartialDerivativeBspline.hpp"
#include "opt/basis/wavelet/boundary/operation/OperationEvalPartialDerivativeWaveletBoundary.hpp"
#include "opt/basis/wavelet/modified/operation/OperationEvalPartialDerivativeModWavelet.hpp"
#include "opt/basis/wavelet/noboundary/operation/OperationEvalPartialDerivativeWavelet.hpp"

namespace sg
{
namespace opt
{

base::OperationEval *createOperationEval(base::Grid &grid)
{
    if (std::strcmp(grid.getType(), "linear") == 0)
    {
        return new OperationEvalLinear(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
    {
        return new OperationEvalLinearBoundary(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)
    {
        return new OperationEvalLinearClenshawCurtis(grid.getStorage(),
                dynamic_cast<LinearClenshawCurtisGrid &>(grid).getCosineTable());
    } else if (std::strcmp(grid.getType(), "modLinear") == 0)
    {
        return new OperationEvalModLinear(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "Bspline") == 0)
    {
        return new OperationEvalBspline(grid.getStorage(),
                dynamic_cast<BsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0)
    {
        return new OperationEvalBsplineBoundary(grid.getStorage(),
                dynamic_cast<BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        return new OperationEvalBsplineClenshawCurtis(grid.getStorage(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getDegree(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getCosineTable());
    } else if (std::strcmp(grid.getType(), "modBspline") == 0)
    {
        return new OperationEvalModBspline(grid.getStorage(),
                dynamic_cast<ModBsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "Wavelet") == 0)
    {
        return new OperationEvalWavelet(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0)
    {
        return new OperationEvalWaveletBoundary(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "modWavelet") == 0) {
        return new OperationEvalModWavelet(grid.getStorage());
    } else
    {
        throw base::factory_exception("OperationEval is not implemented for this grid type.");
    }
}

OperationEvalGradient *createOperationEvalGradient(base::Grid &grid)
{
    if (std::strcmp(grid.getType(), "Bspline") == 0)
    {
        return new OperationEvalGradientBspline(grid.getStorage(),
                dynamic_cast<BsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0)
    {
        return new OperationEvalGradientBsplineBoundary(grid.getStorage(),
                dynamic_cast<BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        return new OperationEvalGradientBsplineClenshawCurtis(grid.getStorage(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getDegree(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getCosineTable());
    } else if (std::strcmp(grid.getType(), "modBspline") == 0)
    {
        return new OperationEvalGradientModBspline(grid.getStorage(),
                dynamic_cast<ModBsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "Wavelet") == 0)
    {
        return new OperationEvalGradientWavelet(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0)
    {
        return new OperationEvalGradientWaveletBoundary(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "modWavelet") == 0)
    {
        return new OperationEvalGradientModWavelet(grid.getStorage());
    } else
    {
        throw base::factory_exception(
                "OperationEvalGradient is not implemented for this grid type.");
    }
}

OperationEvalHessian *createOperationEvalHessian(base::Grid &grid)
{
    if (std::strcmp(grid.getType(), "Bspline") == 0)
    {
        return new OperationEvalHessianBspline(grid.getStorage(),
                dynamic_cast<BsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0)
    {
        return new OperationEvalHessianBsplineBoundary(grid.getStorage(),
                dynamic_cast<BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        return new OperationEvalHessianBsplineClenshawCurtis(grid.getStorage(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getDegree(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getCosineTable());
    } else if (std::strcmp(grid.getType(), "modBspline") == 0)
    {
        return new OperationEvalHessianModBspline(grid.getStorage(),
                dynamic_cast<ModBsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "Wavelet") == 0)
    {
        return new OperationEvalHessianWavelet(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0)
    {
        return new OperationEvalHessianWaveletBoundary(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "modWavelet") == 0)
    {
        return new OperationEvalHessianModWavelet(grid.getStorage());
    } else
    {
        throw base::factory_exception(
                "OperationEvalHessian is not implemented for this grid type.");
    }
}

OperationMultipleHierarchisation *createOperationMultipleHierarchisation(base::Grid &grid)
{
    if (std::strcmp(grid.getType(), "linear") == 0)
    {
        return new OperationMultipleHierarchisationLinear(
                dynamic_cast<LinearGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
    {
        return new OperationMultipleHierarchisationLinearBoundary(
                dynamic_cast<LinearTrapezoidBoundaryGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)
    {
        return new OperationMultipleHierarchisationLinearClenshawCurtis(
                dynamic_cast<LinearClenshawCurtisGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "modLinear") == 0)
    {
        return new OperationMultipleHierarchisationModLinear(
                dynamic_cast<ModLinearGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "Bspline") == 0)
    {
        return new OperationMultipleHierarchisationBspline(
                dynamic_cast<BsplineGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0)
    {
        return new OperationMultipleHierarchisationBsplineBoundary(
                dynamic_cast<BsplineTrapezoidBoundaryGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        return new OperationMultipleHierarchisationBsplineClenshawCurtis(
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "modBspline") == 0)
    {
        return new OperationMultipleHierarchisationModBspline(
                dynamic_cast<ModBsplineGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "Wavelet") == 0)
    {
        return new OperationMultipleHierarchisationWavelet(
                dynamic_cast<WaveletGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0)
    {
        return new OperationMultipleHierarchisationWaveletBoundary(
                dynamic_cast<WaveletTrapezoidBoundaryGrid &>(grid));
    } else if (std::strcmp(grid.getType(), "modWavelet") == 0) {
        return new OperationMultipleHierarchisationModWavelet(
                dynamic_cast<ModWaveletGrid &>(grid));
    } else
    {
        throw base::factory_exception(
                "OperationMultipleHierarchisation is not implemented for this grid type.");
    }
}

OperationEvalPartialDerivative *createOperationEvalPartialDerivative(base::Grid &grid)
{
    if (std::strcmp(grid.getType(), "Bspline") == 0)
    {
        return new OperationEvalPartialDerivativeBspline(grid.getStorage(),
                dynamic_cast<BsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0)
    {
        return new OperationEvalPartialDerivativeBsplineBoundary(grid.getStorage(),
                dynamic_cast<BsplineTrapezoidBoundaryGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0)
    {
        return new OperationEvalPartialDerivativeBsplineClenshawCurtis(grid.getStorage(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getDegree(),
                dynamic_cast<BsplineClenshawCurtisGrid &>(grid).getCosineTable());
    } else if (std::strcmp(grid.getType(), "modBspline") == 0)
    {
        return new OperationEvalPartialDerivativeModBspline(grid.getStorage(),
                dynamic_cast<ModBsplineGrid &>(grid).getDegree());
    } else if (std::strcmp(grid.getType(), "Wavelet") == 0)
    {
        return new OperationEvalPartialDerivativeWavelet(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0)
    {
        return new OperationEvalPartialDerivativeWaveletBoundary(grid.getStorage());
    } else if (std::strcmp(grid.getType(), "modWavelet") == 0) {
        return new OperationEvalPartialDerivativeModWavelet(grid.getStorage());
    } else
    {
        throw base::factory_exception(
                "OperationEvalPartialDerivative is not implemented for this grid type.");
    }
}

}
}
