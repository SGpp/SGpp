// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "PredictiveRefinementIndicator.hpp"
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <map>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {


    PredictiveRefinementIndicator::PredictiveRefinementIndicator(Grid* grid, DataMatrix* dataSet, DataVector* errorVector,
        size_t refinements_num, float_t threshold)

    {
      //find out what type of grid is used;
      gridType = grid->getType();

      //set global Variables accordingly
      this->dataSet = dataSet;
      this->errorVector = errorVector;
      this->refinementsNum = refinements_num;
      this->threshold = threshold;
    }

    float_t PredictiveRefinementIndicator::operator ()(AbstractRefinement::index_type* gridPoint) {
      //calculate the floor and ceiling of the support on dimension of the grid point.
      DataVector floorMask(dataSet->getNcols());
      DataVector ceilingMask(dataSet->getNcols());
      buildGPSupportMask(gridPoint, &floorMask, &ceilingMask);

      //level, index and evaulation of a gridPoint in dimension d
      AbstractRefinement::level_t level = 0;
      AbstractRefinement::index_t index = 0;
      float_t valueInDim;

      //the actuall value of the errorIndicator
      float_t errorIndicator = 0;


      //counter of contributions - for DEBUG purposes
      //size_t counter = 0;

      //go through the whole dataset. -> if data point on the support of the grid point in all dim then calculate error Indicator.
      for (size_t row = 0; row < dataSet->getNrows(); ++row) {
        //if on the support of the grid point in all dim
        if (isOnSupport(&floorMask, &ceilingMask, row)) {
          //counter for DEBUG
          //++counter;
          //calculate error Indicator
          for (size_t dim = 0; dim < gridPoint->dim(); ++dim) {

            level = gridPoint->getLevel(dim);
            index = gridPoint->getIndex(dim);

            valueInDim = dataSet->get(row, dim);

            errorIndicator += ((RefinementFunctor::value_type) basisFunctionEvalHelper(level, index, valueInDim) * errorVector->get(row));
          }
        }
      }

      //DEBUG
      //std::cout << gridPoint->toString() << " with error estimate " << errorIndicator << ",caused by " << counter << "contribs - in average: " << (errorIndicator/static_cast<float_t>(counter)) << "\n";
      return errorIndicator;

    }

    float_t PredictiveRefinementIndicator::operator ()(GridStorage* storage, size_t seq) {
      return errorVector->get(seq);
    }


    float_t PredictiveRefinementIndicator::basisFunctionEvalHelper(AbstractRefinement::level_t level, AbstractRefinement::index_t index, float_t value) {
      if (gridType == base::GridType::Linear) {
        // linear basis
        LinearBasis<AbstractRefinement::level_t, AbstractRefinement::index_t> linBasis;
        return linBasis.eval(level, index, value);
      } else if (gridType == base::GridType::LinearL0Boundary) {
        // linear Basis with Boundaries
        LinearBoundaryBasis<AbstractRefinement::level_t, AbstractRefinement::index_t> linBoundBasis;
        return linBoundBasis.eval(level, index, value);
      } else if (gridType == base::GridType::ModLinear) {
        // modified linear basis
        LinearModifiedBasis<AbstractRefinement::level_t, AbstractRefinement::index_t> modLinBasis;
        return modLinBasis.eval(level, index, value);
      } else {
        // not found.
        return 0.0f;
      }
    }

    size_t PredictiveRefinementIndicator::getRefinementsNum() {
      return refinementsNum;
    }

    float_t PredictiveRefinementIndicator::getRefinementThreshold() {
      return threshold;
    }

    float_t PredictiveRefinementIndicator::start() {
      return 0.0;
    }

    bool PredictiveRefinementIndicator::isOnSupport(
      DataVector* floorMask, DataVector* ceilingMask, size_t row) {

      //go through all cols of the dataset
      //=> go through all samples in dataset and check if in dim "col" "valueInDim" is on support
      for (size_t col = 0; col < dataSet->getNcols(); ++col) {
        float_t valueInDim = dataSet->get(row, col);

        if (valueInDim < floorMask->get(col) || valueInDim >= ceilingMask->get(col) ) {
          return false;
        }
      }

      DataVector vector(dataSet->getNcols());
      dataSet->getRow(row, vector);

      //Debug
      //  std::cout << vector.toString() << " is on support of " << floorMask->toString() << " & " << ceilingMask->toString() << std::endl;

      return true;
    }

    void PredictiveRefinementIndicator::buildGPSupportMask(
      AbstractRefinement::index_type* gridPoint, DataVector* floorMask, DataVector* ceilingMask) {

      AbstractRefinement::level_t level;
      AbstractRefinement::index_t index;

      //in each dimension, get level and index, calculate min and max of supp(GridPointBasisFunction)
      for (size_t dim = 0; dim < gridPoint->dim(); ++dim) {
        level = gridPoint->getLevel(dim);
        index = gridPoint->getIndex(dim);

        floorMask->set(dim, (index - 1.0) / (1 << (level)));
        ceilingMask->set(dim, (index + 1.0) / (1 << (level)));

        //DEBUG
        //    std::cout << "floor: " << floorMask->get(dim) << "ceiling " << ceilingMask->get(dim) <<std::endl;
      }
    }




  } /* namespace base */
} /* namespace SGPP */
