/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef X86SIMDLINEARMULT_H
#define X86SIMDLINEARMULT_H

#include "base/grid/GridStorage.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/ChunkSizes.h"


namespace sg {
namespace parallel {

class X86SimdLinearMult
{
public:
	X86SimdLinearMult(sg::base::DataMatrix* level, sg::base::DataMatrix* index, sg::base::DataMatrix* dataset, sg::base::DataVector& alpha, sg::base::DataVector& result);

	void operator()(size_t start_index_data, size_t end_index_data);
	const size_t chunkGridPoints;
	const size_t chunkDataPoints;

private:
	sg::base::DataMatrix *_level;
	sg::base::DataMatrix *_index;
	sg::base::DataMatrix *_dataset;
	sg::base::DataVector &_alpha;
	sg::base::DataVector &_result;
};

}
}

#endif // X86SIMDLINEARMULT_H
