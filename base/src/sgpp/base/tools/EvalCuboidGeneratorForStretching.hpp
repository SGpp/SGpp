/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef EVALCUBOIDGENERATORFORSTRETCHING_HPP
#define EVALCUBOIDGENERATORFORSTRETCHING_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
//#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>
#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class builds a cuboid in the d-dimension space. This
     * cuboid is discretesized by a regular full grid.
     *
     * @version $HEAD$
     */
    class EvalCuboidGeneratorForStretching {
      private:
        /// number of dimensions
        size_t numDimensions;

        /**
         * This function is a recursive implementation in order the build the evaluation cuboid
         *
         * @param evalPoints vector of dynamic size into which the points are "submitted" during calculation
         * @param curPoint a current point in the d-dimensional space which which is adjusted during this recursive calculations
         * @param myStretching streching
         * @param points number of points used in every dimension
         * @param curDim current dimension in recursive cuboid construction
         */
        void getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, Stretching& myStretching, size_t points, size_t curDim);

      public:
        /**
         * Constructor
         */
        EvalCuboidGeneratorForStretching();

        /**
         * Destructor
         */
        ~EvalCuboidGeneratorForStretching();

        /**
         * This function builds an cuboid which will be stored into the EvaluationPoint
         * variable of this function.
         * This is by done by building a cuboid around a bounding box. Be aware that this
         * function returns points to the power of d points.
         *
         * @param EvaluationPoints DataMatrix that will contain the evaluation points afterwards
         * @param SubDomain the bounding box of the evaluation cuboid
         * @param points number of points used in every dimension
         */
        void getEvaluationCuboid(DataMatrix& EvaluationPoints, Stretching& SubDomain, size_t points);
    };

  }
}

#endif /* EVALCUBOIDGENERATORFORSTRETCHING_HPP */
