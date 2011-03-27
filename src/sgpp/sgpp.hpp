/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGPP_HPP_
#define SGPP_HPP_

#include "grid/GridStorage.hpp"
#include "grid/GridDataBase.hpp"

#include "algorithm/datadriven/AlgorithmDGEMV.hpp"
#include "algorithm/common/GetAffectedBasisFunctions.hpp"

#include "algorithm/datadriven/test_dataset.hpp"
#include "algorithm/datadriven/DMSystemMatrix.hpp"

#include "basis/basis.hpp"

// @todo (heinecke) check if this can be removed
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

#include "grid/Grid.hpp"
#include "grid/common/BoundingBox.hpp"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/generation/StandardGridGenerator.hpp"
#include "grid/generation/BoundaryGridGenerator.hpp"
#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "grid/generation/GridGenerator.hpp"
#include "grid/generation/hashmap/HashGenerator.hpp"
#include "grid/generation/hashmap/HashRefinement.hpp"
#include "grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"


namespace sg
{

typedef linearboundaryBase<unsigned int, unsigned int> SLinearBoundaryBase;
typedef linear_base<unsigned int, unsigned int> SLinearBase;
typedef modified_linear_base<unsigned int, unsigned int> SModLinearBase;

typedef AlgorithmDGEMV<SLinearBase> SGridOperationB;
typedef AlgorithmDGEMV<SLinearBoundaryBase> SGridBoundaryOperationB;
typedef AlgorithmDGEMV<SModLinearBase> SGridModOperationB;

}

#endif /*SGPP_HPP_*/
