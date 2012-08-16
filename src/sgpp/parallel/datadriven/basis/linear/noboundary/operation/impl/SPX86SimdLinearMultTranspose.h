/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPX86SIMDLINEARMULTTRANSPOSE_H
#define SPX86SIMDLINEARMULTTRANSPOSE_H

#include "base/grid/GridStorage.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/ChunkSizes.h"

namespace sg {
namespace parallel {

class SPX86SimdLinearMultTranspose
{
public:
	SPX86SimdLinearMultTranspose(sg::base::DataMatrixSP* level, sg::base::DataMatrixSP* index, sg::base::DataMatrixSP* dataset, sg::base::DataVectorSP& source, sg::base::DataVectorSP& result);

	void operator()(size_t start_index_grid, size_t end_index_grid);
	const size_t chunkGridPoints;
	const size_t chunkDataPoints;

private:
	sg::base::DataMatrixSP *_level;
	sg::base::DataMatrixSP *_index;
	sg::base::DataMatrixSP *_dataset;
	sg::base::DataVectorSP &_source;
	sg::base::DataVectorSP &_result;
};

}
}

#endif // SPX86SIMDLINEARMULTTRANSPOSE_H
