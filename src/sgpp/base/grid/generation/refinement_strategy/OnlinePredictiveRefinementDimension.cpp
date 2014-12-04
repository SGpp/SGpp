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

#include "OnlinePredictiveRefinementDimension.hpp"
#include "base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp"
#include "base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp"


namespace sg
{
namespace base
{

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

    Grid* predictiveGrid = Grid::createLinearGrid(dim);
    GridGenerator* predictiveGridGenerator = predictiveGrid->createGridGenerator();
    predictiveGridGenerator->regular(2);

    GridStorage* predictiveGridStorage = predictiveGrid->getStorage();
    size_t predictiveGridSize = predictiveGrid->getSize();

    OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*predictiveGrid, trainDataset);

    // Display points on predictive Grid

    GridStorage::grid_map_iterator predictive_end_iter = predictiveGrid->getStorage()->end();
    for (GridStorage::grid_map_iterator iter = predictiveGrid->getStorage()->begin();
            iter != predictive_end_iter; iter++)
    {
        // std::cout << predictiveGrid->getStorage()->seq(iter->first) << ": " << iter->first->toString() << std::endl;
    }


    GridIndex* gridIndex;
    GridStorage::grid_map_iterator end_iter = storage->end();
    for (GridStorage::grid_map_iterator iter = storage->begin();
            iter != end_iter; iter++)
    {

        gridIndex = iter->first;

        //std::cout << "Point: " << storage->seq(gridIndex) << " (" << gridIndex->toString() << ")" << std::endl;

        // Refinability

        bool refinable = false;
        for( size_t k=0; k <dim; k++)
        {
            if( !hasLeftChild(storage, gridIndex, k))
            {
                refinable = true;
            }
        }
        if (!refinable)
        {
            continue;
        }

        // Bounding box = support of refinable point

        // std::cout << "New bounding box:" << std::endl;

        DimensionBoundary* boundaries = new DimensionBoundary[dim];

        for (size_t k = 0; k < dim; k++ )
        {
            // std::cout << "Dimension " << k << std::endl;

            unsigned int index = gridIndex->getIndex(k);
            unsigned int level = gridIndex->getLevel(k);
            double intval = pow(2.0, -static_cast<double>(level));

            DimensionBoundary boundary;
            boundary.leftBoundary = (index-1) * intval;
            boundary.rightBoundary = (index+1) * intval;
            boundaries[k] = boundary;

            // std::cout << "left boundary: " << boundary.leftBoundary << std::endl;
            // std::cout << "right boundary: " << boundary.rightBoundary << std::endl;
        }

        BoundingBox* bb = new BoundingBox(dim, boundaries);
        predictiveGrid->setBoundingBox(*bb);

        // All numerators
        DataVector numerators(predictiveGridSize);
        eval->multTranspose(*errors, numerators);
        //std::cout << "Numerators: " << numerators.toString() << std::endl;
        //std::cout << "Errors: " << errors->toString() << std::endl;
        numerators.sqr();

        // All denominators
        DataVector denominators(predictiveGridSize);

        for (size_t j = 0; j < predictiveGridSize; j++)
        {
            DataVector single(predictiveGridSize);
            single.setAll(0.0);
            single.set(j, 1.0);

            DataVector col(numData);
            eval->mult(single, col);

            col.sqr();
            denominators.set(j, col.sum());
        }

        // Calculate the indicator value for all refinable dimensions

        for (size_t k = 0; k < dim; k++ )
        {
            // std::cout << "Dimension: " << k << std::endl;

            // Is the current point refinable in the dimension k?
            if( hasLeftChild(storage, gridIndex, k) || hasRightChild(storage, gridIndex, k) )
            {
                continue;
            }

            // Identify the index of the two grid points
            // and calculate their indicator value

            GridIndex childIndex(dim);
            size_t childSeq;
            double value1 = 0;
            double value2 = 0;

            // Left Child
            childIndex.set(k, 2, 1);

            // other dimensions
            for(size_t j=0;j<dim;j++)
            {
                if (j != k)
                {
                    childIndex.set(j, 1, 1);
                }
            }
            childSeq = predictiveGridStorage->seq(&childIndex);

            if (denominators.get(childSeq) != 0)
            {
                value1 = numerators.get(childSeq) / denominators.get(childSeq);
                //std::cout << "Left Child: sequence: " << childSeq << std::endl;
                //std::cout << "Left Child: numerator: " << numerators.get(childSeq) << std::endl;
                //std::cout << "Left Child: denominators: " << denominators.get(childSeq) << std::endl;
            }
            else
            {
                // No points or only boundary points
                //std::cout << "Left Child: has no points" << std::endl;
                value1 = 0;
            }

            /*
            std::cout << "Left Child: " << value1 << std::endl;
            std::cout << numerators.get(childSeq) << std::endl;
            std::cout << denominators.get(childSeq) << std::endl;
			*/

            // Right Child
            childIndex.set(k, 2, 3);
            for(size_t j=0;j<dim;j++)
            {
                if (j != k)
                {
                    childIndex.set(j, 1, 1);
                }
            }
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

            /*
            std::cout << "Right Child: " << value2 << std::endl;
            std::cout << numerators.get(childSeq) << std::endl;
            std::cout << denominators.get(childSeq) << std::endl;
			*/

            result->insert(std::pair<refinement_key, double>(
                               refinement_key(storage->seq(gridIndex),
                                              (unsigned int) k), value1 + value2)
                          );

            //std::cout << "Point " << storage->seq(gridIndex) << " (dim: " << k << "): " << (value1 + value2) << std::endl;
        }

        delete bb;
        delete[] boundaries;

    }

    // Output result
    /*std::cout << "Result1: " << std::endl;
    for(refinement_map::iterator it = result->begin(); it != result->end(); it++ )
    {
    	std::cout << "(" << it->first.first << ", " << it->first.second << "): " << it->second << std::endl;
    }*/

    if(refinements_num != 0)
    {
        typedef std::vector<std::pair<double, refinement_key> > TmpVec;
        TmpVec v;

        for(refinement_map::iterator it = result->begin(); it != result->end(); it++ )
        {
            v.push_back(std::make_pair(it->second, it->first));
        }

        std::sort(v.begin(), v.end(), std::greater<std::pair<double, refinement_key> >());

        refinement_map result2 = refinement_map();

        size_t i = 0;
        for(TmpVec::iterator it = v.begin(); it != v.end(); it++)
        {
            result2[it->second] = it->first;
            i++;
            if( i == refinements_num )
            {
                break;
            }
        }

        *result = result2;

        // Output result2
        /*
        std::cout << "Result2:" << std::endl;
        for(refinement_map::iterator it = result2.begin(); it != result2.end(); it++ )
        {
        	std::cout << "(" << it->first.first << ", " << it->first.second << "): " << it->second << std::endl;
        }*/
    }

    delete predictiveGrid;
    delete predictiveGridGenerator;
    delete eval;

}

bool refinementPairCompare(const std::pair<OnlinePredictiveRefinementDimension::key_type, OnlinePredictiveRefinementDimension::value_type>& firstEl,
                           const std::pair<OnlinePredictiveRefinementDimension::key_type, OnlinePredictiveRefinementDimension::value_type>& secondEl)
{
    return firstEl.second > secondEl.second;
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

    	if ( val < functor->getRefinementThreshold() ) {
    		continue;
    	}

    	index_type index = storage->get(seq);
    	index.setLeaf(false);

    	// if does not work, try
    	// (*storage)[max_index]->setLeaf(false);

    	printf("Refine: (%d,%d) with value %2.10f\n", (int)seq, (int)dim, val);
    	fflush(stdout);

    	this->refineGridpoint1D(storage, index, dim);
    }
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
