// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "PredictiveRefinementDimensionIndicator.hpp"
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <map>
#include <string>
#include <utility>
#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {


    PredictiveRefinementDimensionIndicator::PredictiveRefinementDimensionIndicator(Grid* grid, DataMatrix* dataSet, DataVector* errorVector,
        size_t refinements_num, float_t threshold, long unsigned int minSupportPoints): minSupportPoints_(minSupportPoints) {
      //find out what type of grid is used;
      gridType = grid->getType();

      //set global Variables accordingly
      this->dataSet = dataSet;
      this->errorVector = errorVector;
      this->refinementsNum = refinements_num;
      this->threshold = threshold;
      this->grid_ = grid;
    }

    float_t PredictiveRefinementDimensionIndicator::operator ()(AbstractRefinement::index_type* gridPoint) {
      //the actuall value of the errorIndicator
      float_t errorIndicator = 0.0;
      float_t denominator = 0.0;
      float_t r22 = 0.0;
      float_t r2phi = 0.0;


      //counter of contributions - for DEBUG purposes
      size_t counter = 0;

      SBasis& basis = const_cast<SBasis&>(grid_->getBasis());
      //go through the whole dataset. -> if data point on the support of the grid point in all dim then calculate error Indicator.
      #pragma omp parallel for schedule(static) reduction(+:errorIndicator,denominator,r22,r2phi,counter)

      for (size_t row = 0; row < dataSet->getNrows(); ++row) {
        //level, index and evaulation of a gridPoint in dimension d
        AbstractRefinement::level_t level = 0;
        AbstractRefinement::index_t index = 0;
        float_t valueInDim;
        //if on the support of the grid point in all dim
        //if(isOnSupport(&floorMask,&ceilingMask,row))
        //{
        //counter for DEBUG
        //++counter;*****
        float_t funcval = 1.0;

        //calculate error Indicator
        for (size_t dim = 0; dim < gridPoint->dim() && funcval != 0; ++dim) {

          level = gridPoint->getLevel(dim);
          index = gridPoint->getIndex(dim);

          valueInDim = dataSet->get(row, dim);

          funcval *=  std::max(float_t(0), basis.eval(level,
                               index,
                               valueInDim));

          //basisFunctionEvalHelper(level,index,valueInDim);
        }

        errorIndicator += funcval * errorVector->get(row)/**errorVector->get(row)*/;
        r22 += errorVector->get(row) * errorVector->get(row);
        r2phi += funcval * errorVector->get(row);
        denominator += funcval * funcval;

        if (funcval != 0.0) counter++;

        //}
      }

      AbstractRefinement::index_type idx(*gridPoint);
      countersMap[idx] = counter;

      if (denominator != 0 && counter >= minSupportPoints_) {
        // to match with OnlineRefDim, use this:
        //return (errorIndicator * errorIndicator) / denominator;

        float_t a = (errorIndicator / denominator);
        return /*r2phi/denominator*/ /*2*r22 - 2*a*r2phi + a*a*denominator*/ a * (2 * r2phi - a * denominator);
        //return fabs(a);
      } else {
        return 0.0;
      }

    }


    float_t PredictiveRefinementDimensionIndicator::operator ()(GridStorage* storage, size_t seq) {
      return errorVector->get(seq);
    }


    float_t PredictiveRefinementDimensionIndicator::runOperator(GridStorage* storage, size_t seq) {
      return (*this)(storage->get(seq));
    }


    float_t PredictiveRefinementDimensionIndicator::basisFunctionEvalHelper(AbstractRefinement::level_t level, AbstractRefinement::index_t index, float_t value) {
      if (gridType == base::GridType::Linear) {
        // linear basis
        LinearBasis<AbstractRefinement::level_t, AbstractRefinement::index_t> linBasis;
        return linBasis.eval(level, index, value);
      } else if (gridType == base::GridType::LinearBoundary ||
          gridType == base::GridType::LinearL0Boundary) {
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

    size_t PredictiveRefinementDimensionIndicator::getRefinementsNum() {
      return refinementsNum;
    }

    float_t PredictiveRefinementDimensionIndicator::getRefinementThreshold() {
      return threshold;
    }

    float_t PredictiveRefinementDimensionIndicator::start() {
      return 0.0;
    }

    bool PredictiveRefinementDimensionIndicator::isOnSupport(
      DataVector* floorMask, DataVector* ceilingMask, size_t row) {

      //go through all cols of the dataset
      //=> go through all samples in dataset and check if in dim "col" "valueInDim" is on support
      for (size_t col = 0; col < dataSet->getNcols(); ++col) {
        float_t valueInDim = dataSet->get(row, col);

        if (valueInDim < floorMask->get(col) || valueInDim >= ceilingMask->get(col) ) {
          return false;
        }
      }

      return true;
    }

    void PredictiveRefinementDimensionIndicator::buildGPSupportMask(
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
