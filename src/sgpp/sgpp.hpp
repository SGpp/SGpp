/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef SGPP_HPP_
#define SGPP_HPP_

// Optimizations currently useless
//#define SGPP_OPTIMIZE

#include "algorithm/classification/AlgorithmDGEMV.hpp"
#include "algorithm/classification/AlgorithmDGEMVBoundaries.hpp"
#include "algorithm/common/GetAffectedBasisFunctions.hpp"
#include "algorithm/common/GetAffectedBasisFunctionsBoundaries.hpp"
#include "algorithm/common/GetAffectedBasisFunctionsBBBoundaries.hpp"
#include "algorithm/classification/test_dataset.hpp"
#include "algorithm/classification/test_dataset_boundary.hpp"
#include "algorithm/classification/DMSystemMatrix.hpp"
#include "algorithm/finance/BlackScholesTimestepMatrix.hpp"

#include "application/finance/BlackScholesSolver.hpp"

#include "basis/basis.hpp"

#include "basis/linear/operation/classification/OperationLaplaceLinear.hpp"
#include "basis/linearboundary/operation/classification/OperationLaplaceLinearBoundary.hpp"
#include "basis/lineartrapezoidboundary/operation/classification/OperationLaplaceLinearTrapezoidBoundary.hpp"
#include "basis/modlinear/operation/classification/OperationLaplaceModLinear.hpp"

#include "data/DataVector.hpp"

#include "grid/Grid.hpp"
#include "grid/GridStorage.hpp"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/generation/StandardGridGenerator.hpp"
#include "grid/generation/BoundaryGridGenerator.hpp"
#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "grid/generation/GridGenerator.hpp"
#include "grid/generation/hashmap/HashGenerator.hpp"
#include "grid/generation/hashmap/HashRefinement.hpp"
#include "grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"

#include "solver/sle/ConjugateGradients.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "solver/ode/ExplicitEuler.hpp"
#include "solver/ode/CrankNicolson.hpp"

#include "tools/finance/IOToolBonnSG.hpp"

#include "linearSystem/LinearSystem.hpp"

namespace sg
{

typedef lineartrapezoidboundaryBase<unsigned int, unsigned int> SLinearTrapezoidBoundaryBase;
typedef linearboundaryBase<unsigned int, unsigned int> SLinearBoundaryBase;
typedef linear_base<unsigned int, unsigned int> SLinearBase;
typedef modified_linear_base<unsigned int, unsigned int> SModLinearBase;
typedef poly_base<unsigned int, unsigned int> SPolyBase;
typedef modified_poly_base<unsigned int, unsigned int> SModPolyBase;

typedef AlgorithmDGEMV<SLinearBase> SGridOperationB;
typedef AlgorithmDGEMV<SLinearBoundaryBase> SGridBoundaryOperationB;
typedef AlgorithmDGEMV<SLinearTrapezoidBoundaryBase> SGridTrapezoidBoundaryOperationB;
typedef AlgorithmDGEMV<SModLinearBase> SGridModOperationB;

}

#endif /*SGPP_HPP_*/
