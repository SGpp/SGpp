/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Maxim Schmidt (maxim.schmidt@tum.de)
#include <algorithm>
#include <limits>
#include <list>
#include <vector>
#include <utility>
#include <cmath>

#include "OnlinePredictiveRefinementDimension.hpp"
#include <sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp>
#include <sgpp/base/basis/linear/noboundary/LinearBasis.hpp>
#include <sgpp/parallel/operation/ParallelOpFactory.hpp>


namespace sg
{
namespace base
{

bool refinementPairCompare(const std::pair<OnlinePredictiveRefinementDimension::key_type, OnlinePredictiveRefinementDimension::value_type>& firstEl,
                           const std::pair<OnlinePredictiveRefinementDimension::key_type, OnlinePredictiveRefinementDimension::value_type>& secondEl)
{
    return firstEl.second > secondEl.second;
}

const sg::parallel::VectorizationType OnlinePredictiveRefinementDimension::vecType_ = sg::parallel::VectorizationType::X86SIMD;


void OnlinePredictiveRefinementDimension::collectRefinablePoints(
    GridStorage* storage, size_t refinements_num, refinement_map* result)
{

    using namespace sg::base;

    if (trainDataset == NULL || errors == NULL)
    {
        throw base::application_exception(
            "Training dataset, classes or errors not set");
    }

    size_t dim = storage->dim();
    size_t numData = trainDataset->getNrows();

    //std::cout << "Dataset: " << trainDataset->toString() << std::endl;

    GridStorage* predictiveGridStorage = predictiveGrid_->getStorage();
    size_t predictiveGridSize = predictiveGrid_->getSize();

    OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*predictiveGrid_, *trainDataset);
//    sg::parallel::OperationMultipleEvalVectorized* eval =
//    		sg::op_factory::createOperationMultipleEvalVectorized(*predictiveGrid_,
//                vecType_, trainDataset);

    GridIndex* gridIndex;
    for (GridStorage::grid_map_iterator iter = storage->begin();
            iter != storage->end(); iter++)
    {

        gridIndex = iter->first;

        // Refinability

        bool refinable = false;
        for( size_t k=0; k <dim; k++)
        {
            if( !hasLeftChild(storage, gridIndex, k))
            {
                refinable = true;
                break;
            }
        }
        if (!refinable)
        {
            continue;
        }

        // Bounding box = support of refinable point
        DimensionBoundary* boundaries = new DimensionBoundary[dim];

		//#pragma omp parallel for schedule(static)
        for (size_t k = 0; k < dim; k++ )
        {
            double index = gridIndex->getIndex(k);
            double level = gridIndex->getLevel(k);
            double intval = pow(2.0, -level);
            //std::cout << "level " << level << " intval " << intval << std::endl;

            DimensionBoundary boundary;
            boundary.leftBoundary = (index-1) * intval;
            boundary.rightBoundary = (index+1) * intval;
            boundaries[k] = boundary;
        }

        BoundingBox* bb = new BoundingBox(dim, boundaries);
        predictiveGrid_->setBoundingBox(*bb);

        //std::cout << predictiveGrid_->getBoundingBox()->toString() << std::endl;


        // All numerators
        // FIXME: should be initiated only once
        DataVector numerators(predictiveGridSize);
        eval->multTranspose(*errors, numerators);
        //numerators.sqr();

        // All denominators
        // FIXME: should be initiated only once
        DataVector denominators(predictiveGridSize);
        denominators.setAll(0.0);

//        DataVector single(predictiveGridSize);
//		single.setAll(0.0);
        SBasis& basis = const_cast<SBasis&>(predictiveGrid_->getBasis());
		#pragma omp parallel for schedule(static)
        for (size_t j = 0; j < predictiveGridSize; j++)
        {

//            single.set(j, 1.0);
//            DataVector  private(single)col(numData);
//            eval->mult(single, col);
//            col.sqr();
//            denominators.set(j, col.sum());
//            single.set(j, 0.0);


            GridIndex predGridIdx = predictiveGridStorage->get(j);
            for (size_t point_idx = 0; point_idx < numData; point_idx++){
            	double prod = 1.0;
            	for (size_t d = 0; d < dim && prod != 0; d++){
                	double transformedCoords = (trainDataset->get(point_idx, d) -
                			boundaries[d].leftBoundary)/(boundaries[d].rightBoundary - boundaries[d].leftBoundary);
            		prod *= std::max(0.0, basis.eval(predGridIdx.getLevel(d),
            				predGridIdx.getIndex(d),
            				transformedCoords));
            	}
            	denominators[j] += prod*prod;

            }


        }


        // Calculate the indicator value for all refinable dimensions
        GridIndex childIndex(dim);
        for(size_t j=0;j<dim;j++)
		{
			childIndex.set(j, 1, 1);
		}
        for (size_t k = 0; k < dim; k++ )
        {
            // Is the current point refinable in the dimension k?
            if( hasLeftChild(storage, gridIndex, k) || hasRightChild(storage, gridIndex, k) )
            {
                continue;
            }

            // Identify the index of the two grid points
            // and calculate their indicator value

            size_t childSeq;

            double value1 = 0;
            double value2 = 0;

            // Left Child
            childIndex.set(k, 2, 1);

            childSeq = predictiveGridStorage->seq(&childIndex);

            if (denominators.get(childSeq) != 0)
            {
                value1 = numerators.get(childSeq) / denominators.get(childSeq);
            }
            else
            {
                // No points or only boundary points
                value1 = 0;
            }

            // Right Child
            childIndex.set(k, 2, 3);
            childSeq = predictiveGridStorage->seq(&childIndex);

            if (denominators.get(childSeq) != 0)
            {
                value2 = numerators.get(childSeq) / denominators.get(childSeq);
            }
            else
            {
                // No points or only boundary points
                value2 = 0;
            }

            result->insert(std::pair<refinement_key, double>(
                               refinement_key(storage->seq(gridIndex),
                                              (unsigned int) k), value1 + value2)
                          );

            childIndex.set(k, 1, 1);
        }

        delete bb;
        delete[] boundaries;

    }

    if(refinements_num != 0)
    {
        typedef std::vector<std::pair<key_type, value_type> > TmpVec;
        TmpVec errorsVector;

        /*for(refinement_map::iterator it = result->begin(); it != result->end(); it++ )
        {
            errorsVector.push_back(std::make_pair(it->first, it->second));
        }*/
        std::copy(result->begin(),
        		result->end(),
        	       std::back_inserter<std::vector<std::pair<key_type, value_type> > >(errorsVector));

        std::nth_element(errorsVector.begin(), errorsVector.begin()+refinements_num,  errorsVector.end(), refinementPairCompare);

        //refinement_map result2 = refinement_map();
        result->clear();

        for(TmpVec::iterator it = errorsVector.begin(); it != errorsVector.begin()+std::min(refinements_num, errorsVector.size()); it++)
        {
            (*result)[it->first] = it->second;
        }
    }

    delete eval;

}


void OnlinePredictiveRefinementDimension::refineGridpointsCollection(
    GridStorage* storage, RefinementFunctor* functor,
    size_t refinements_num, size_t* max_indices,
    PredictiveRefinementDimensionIndicator::value_type* max_values)
{}

void OnlinePredictiveRefinementDimension::free_refine(GridStorage* storage,
        RefinementFunctor* functor)
{
    if (storage->size() == 0)
    {
        throw generation_exception("storage empty");
    }

    size_t refinements_num = functor->getRefinementsNum();

    refinement_map result = refinement_map();

    OnlinePredictiveRefinementDimension::collectRefinablePoints(storage, refinements_num, &result);

    for(refinement_map::iterator it = result.begin(); it != result.end(); it++)
    {
    	size_t seq = it->first.first;
    	unsigned int dim = it->first.second;
    	double val = it->second;

    	if ( val < functor->getRefinementThreshold() || val <= functor->start() ) {
    		continue;
    	}

    	(*storage)[it->first.first]->setLeaf(false);

    	printf("Refine: (%d,%d) with value %2.10f\n", (int)seq, (int)dim, val);
    	fflush(stdout);

    	this->refineGridpoint1D(storage, *(*storage)[it->first.first], dim);
    }
    storage->recalcLeafProperty();
}

bool OnlinePredictiveRefinementDimension::hasLeftChild(GridStorage* storage, GridIndex* gridIndex, size_t dim)
{
    GridIndex tmp = GridIndex(*gridIndex);
    storage->left_child(&tmp, dim);
    return storage->has_key(&tmp);
}

bool OnlinePredictiveRefinementDimension::hasRightChild(GridStorage* storage, GridIndex* gridIndex, size_t dim)
{
    GridIndex tmp = GridIndex(*gridIndex);
    storage->right_child(&tmp, dim);
    return storage->has_key(&tmp);
}

size_t OnlinePredictiveRefinementDimension::getIndexOfMin(
    PredictiveRefinementDimensionIndicator::value_type* array,
    size_t length)
{
    size_t min_idx = 0;

    for (size_t i = 1; i < length; i++)
    {
        if (array[i].second < array[min_idx].second)
            min_idx = i;
    }

    return min_idx;
}

void OnlinePredictiveRefinementDimension::setTrainDataset(
    DataMatrix* trainDataset_)
{
    trainDataset = trainDataset_;
}

void OnlinePredictiveRefinementDimension::setErrors(DataVector* errors_)
{
    errors = errors_;
}

double OnlinePredictiveRefinementDimension::basisFunctionEvalHelper(unsigned int level, unsigned int index, double value)
{
    LinearBasis<unsigned int,unsigned int> linBasis;
    return linBasis.eval(level,index,value);
}


} /* namespace base */
} /* namespace sg */
