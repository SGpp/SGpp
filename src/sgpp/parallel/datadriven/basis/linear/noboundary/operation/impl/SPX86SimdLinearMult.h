/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPX86SIMDLINEARMULT_H
#define SPX86SIMDLINEARMULT_H

#include "base/grid/GridStorage.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/ChunkSizes.h"

namespace sg {
namespace parallel {

class SPX86SimdLinearMult
{
public:
	SPX86SimdLinearMult(sg::base::DataMatrixSP* level, sg::base::DataMatrixSP* index, sg::base::DataMatrixSP* dataset, sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result);

	void operator()(size_t start_index_data, size_t end_index_data);
	const size_t chunkGridPoints;
	const size_t chunkDataPoints;

private:
	sg::base::DataMatrixSP *_level;
	sg::base::DataMatrixSP *_index;
	sg::base::DataMatrixSP *_dataset;
	sg::base::DataVectorSP &_alpha;
	sg::base::DataVectorSP &_result;
};

}
}
#endif // SPX86SIMDLINEARMULT_H
