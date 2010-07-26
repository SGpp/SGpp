/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ALGORITHMEVALUATIONITERATIVE_HPP
#define ALGORITHMEVALUATIONITERATIVE_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"
#include "basis/modwavelet/modified_wavelet_base.hpp"
#include "basis/modbspline/modified_bspline_base.hpp"
#include "basis/linear/boundary/linearboundaryBase.hpp"

#include <vector>
#include <utility>

namespace sg {

/**
 * Iterative algorithms for sparse grid
 * function evaluation
 */
template<class BASIS>
class AlgorithmEvaluationIterative {
public:
	AlgorithmEvaluationIterative(GridStorage* storage) :
        storage(storage) {
    }

    ~AlgorithmEvaluationIterative() {
    }

    double operator()(BASIS& basis, std::vector<double>& point, DataVector& alpha) {
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;
        typedef GridStorage::index_pointer index;

        double result = 0.0;
        double curSupport;
        size_t dims = this->storage->dim();
        size_t storageSize = this->storage->size();
        level_type curLevel;
        index_type curIndex;
        index curGridPoint;

        for (size_t i = 0; i < storageSize; i++)
        {
        	curSupport = 1.0;
        	curGridPoint = (*this->storage)[i];

        	for (size_t d = 0; d < dims; d++)
        	{
        		curGridPoint->get(d, curLevel, curIndex);
        		curSupport *= basis.evalSave(curLevel, curIndex, point[d]);
        	}

        	result += (curSupport * alpha[i]);
        }

        return result;
    }

protected:
    GridStorage* storage;
};

}

#endif /* ALGORITHMEVALUATIONITERATIVE_HPP */
