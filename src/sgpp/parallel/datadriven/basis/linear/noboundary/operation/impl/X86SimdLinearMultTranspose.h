/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef X86SIMDLINEARMULTTRANSPOSE_H
#define X86SIMDLINEARMULTTRANSPOSE_H

#include "base/grid/GridStorage.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/ChunkSizes.h"

namespace sg {
namespace parallel {

class X86SimdLinearMultTranspose
{
public:
	X86SimdLinearMultTranspose(sg::base::DataMatrix* level, sg::base::DataMatrix* index, sg::base::DataMatrix* dataset, sg::base::DataVector& source, sg::base::DataVector& result);

	void operator()(size_t start_index_grid, size_t end_index_grid);
	const size_t chunkGridPoints;
	const size_t chunkDataPoints;

public:
	sg::base::DataMatrix *_level;
	sg::base::DataMatrix *_index;
	sg::base::DataMatrix *_dataset;
	sg::base::DataVector &_source;
	sg::base::DataVector &_result;
};

}
}

#endif // X86SIMDMULTTRANSPOSE_H
