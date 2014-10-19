/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
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
    	std::cout << predictiveGrid->getStorage()->seq(iter->first) << ": " << iter->first->toString() << std::endl;
    }


    GridIndex* gridIndex;
    GridStorage::grid_map_iterator end_iter = storage->end();
    for (GridStorage::grid_map_iterator iter = storage->begin();
            iter != end_iter; iter++)
    {

        gridIndex = iter->first;

        std::cout << "Point: " << storage->seq(gridIndex) << " (" << gridIndex->toString() << ")" << std::endl;

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

            //std::cout << "left boundary: " << boundary.leftBoundary << std::endl;
            //std::cout << "right boundary: " << boundary.rightBoundary << std::endl;
        }

        BoundingBox* bb = new BoundingBox(dim, boundaries);
        predictiveGrid->setBoundingBox(*bb);

        // All numerators
        DataVector numerators(predictiveGridSize);
        eval->multTranspose(*errors, numerators);
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
        	std::cout << "Dimension: " << k << std::endl;

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
            for(size_t j=0;j<dim;j++)
            {
                if (j != k)
                {
                    childIndex.set(j, 1, 1);
                }
            }
            childSeq = predictiveGridStorage->seq(&childIndex);

            if (denominators.get(childSeq) != 0) {
            	value1 = numerators.get(childSeq) / denominators.get(childSeq);
            } else {
            	// No points or only boundary points
            	value1 = 0;
            }

           // std::cout << numerators.get(childSeq) << std::endl;
            //std::cout << denominators.get(childSeq) << std::endl;
            std::cout << "Left Child: " << value1 << std::endl;

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

            if (denominators.get(childSeq) != 0) {
            	value2 = numerators.get(childSeq) / denominators.get(childSeq);
            } else {
            	// No points or only boundary points
            	value2 = 0;
            }

            // TODO only the highest
            result->insert(std::pair<refinement_key, double>(
                               refinement_key(storage->seq(gridIndex),
                                              (unsigned int) k), value1 + value2)
                          );

           // std::cout << numerators.get(childSeq) << std::endl;
            //std::cout << denominators.get(childSeq) << std::endl;
            std::cout << "Right Child: " << value2 << std::endl;

            std::cout << std::endl;
        }

        delete bb;
        delete[] boundaries;

    }

    delete predictiveGrid;
    delete predictiveGridGenerator;
    delete eval;


    /* ======================================== */

    /*

    //this refinement algorithm uses the predictive refinement indicator.
    //dynamic casting is used to maintain the signature of the algorithm, but still be able to use the
    //predictive refinement indicator with it.
    PredictiveRefinementDimensionIndicator* errorIndicator =
    		dynamic_cast<PredictiveRefinementDimensionIndicator*>(functor);

    //size_t min_idx = 0;
    std::vector<value_type> errors;

    // max value equals min value
    //PredictiveRefinementDimensionIndicator::value_type* max_values =static_cast<PredictiveRefinementDimensionIndicator::value_type*>(mv);
    //PredictiveRefinementDimensionIndicator::value_type max_value = max_values[min_idx];

    index_type index;
    GridStorage::grid_map_iterator end_iter = storage->end();


    // start iterating over whole grid
    for (GridStorage::grid_map_iterator iter = storage->begin();
    		iter != end_iter; iter++) {
    	index = *(iter->first);

    	GridStorage::grid_map_iterator child_iter;

    	// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
    	// if yes, check whether it belongs to the refinements_num largest ones
    	for (size_t d = 0; d < storage->dim(); d++) {
    		index_t source_index;
    		level_t source_level;
    		index.get(d, source_level, source_index);
    		double error = functor->start();
    		//errorIndicator->setActiveDim(d);

    		// test existence of left child
    		index.set(d, source_level + 1, 2 * source_index - 1);
    		child_iter = storage->find(&index);

    		if (child_iter == end_iter) {
    			//use the predictive error indicator, which takes a pointer to the grid point object
    			//instead of the storage index
    			error += (*errorIndicator)(&index);
    		}

    		// test existance of right child
    		index.set(d, source_level + 1, 2 * source_index + 1);
    		child_iter = storage->find(&index);

    		if (child_iter == end_iter) {//use predictive refinement indicator
    			//use the predictive error indicator, which takes a pointer to the grid point object
    			//instead of the storage index
    			error += (*errorIndicator)(&index);
    		}
    		// reset current grid point in dimension d
    		index.set(d, source_level, source_index);

    		if (error > functor->start()){
    			errors.push_back(error);
    		}

    		if (error > iThreshold_) {
    			key_type key(storage->seq(iter->first), d);
    //				if (refinementCollection_.find(key)
    //						!= refinementCollection_.end()) {
    //					refinementCollection_[key] = refinementCollection_[key] + error;
    //				} else {
    				refinementCollection_[key] = error;
    //				}
    		}
    	}
}

    // get new thershold for the next iteration
    std::nth_element(errors.begin(), errors.begin()+refinements_num, errors.end(), doubleReverseCompare);
    //iThreshold_ = errors[refinements_num];
     *
     */

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
{
    /*
       PredictiveRefinementDimensionIndicator::value_type max_value;
       //size_t max_index;
       //PredictiveRefinementDimensionIndicator::value_type* max_values =static_cast<PredictiveRefinementDimensionIndicator::value_type*>(mv);
       // now refine all grid points which satisfy the refinement criteria
       double threshold = functor->getRefinementThreshold();

       std::vector< std::pair<key_type, value_type> > errorsVector;
       std::cout << "refinement collection size " << refinementCollection_.size() << std::endl;
       std::copy(refinementCollection_.begin(),
                 refinementCollection_.end(),
                 std::back_inserter<std::vector<std::pair<key_type, value_type> > >(errorsVector));

       std::nth_element(errorsVector.begin(), errorsVector.begin()+refinements_num,
                        errorsVector.end(), refinementPairCompare);

       std::vector<std::pair<key_type, value_type> >::const_iterator iter;
       for (iter = errorsVector.begin(); iter < errorsVector.begin() + refinements_num; iter++)
       {
           if (iter->second > functor->start() && iter->second >= threshold)
           {
               index_type index((*storage)[iter->first.first]);
               //Sets leaf property of index, which is refined to false
               (*storage)[iter->first.first]->setLeaf(false);
               std::cout << "Refining grid point " << iter->first.first << " dim " << iter->first.second << " value "
               << iter->second << std::endl;
               this->refineGridpoint1D(storage, index, iter->first.second);
           }
       }
    */
    /*for (size_t i = 0; i < refinements_num; i++) {
    	max_value = max_values[i];
    	max_index = max_indices[i];

    	if (max_value.second > functor->start()
    			&& max_value.second >= threshold) {
    		index_type index((*storage)[max_index]);
    		//Sets leaf property of index, which is refined to false
    		(*storage)[max_index]->setLeaf(false);
    		std::cout << "Refining grid point " << max_index << " value "
    				<< max_value.second << std::endl;
    		this->refineGridpoint1D(storage, index, max_value.first);
    	}
}*/

}

void OnlinePredictiveRefinementDimension::free_refine(GridStorage* storage,
        PredictiveRefinementDimensionIndicator* functor)
{
    /*
       if (storage->size() == 0)
       {
           throw generation_exception("storage empty");
       }

       // the functor->getRefinementsNum() largest grid points should be refined.
       // gather them in an array max_values
       size_t refinements_num = functor->getRefinementsNum();
       // values
       PredictiveRefinementDimensionIndicator::value_type* max_values =
           new PredictiveRefinementDimensionIndicator::value_type[refinements_num];
       // indices
       size_t* max_indices = new size_t[refinements_num];

       // initialization
       for (size_t i = 0; i < refinements_num; i++)
       {
           max_values[i].second = functor->start();
           max_indices[i] = 0;
       }

       //collectRefinablePoints(storage, functor, refinements_num, max_indices,
       //                       max_values);
       // now refine all grid points which satisfy the refinement criteria
       refineGridpointsCollection(storage, functor, refinements_num, max_indices,
                                  max_values);
       refinementCollection_.clear();
       delete[] max_values;
       delete[] max_indices;
    */
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
