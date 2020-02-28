// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/OperationEvalModLinear.hpp>



#include <sgpp/globaldef.hpp>

#include <utility>
#include <vector>


namespace sgpp {
namespace base {

double OperationEvalModLinear::eval(const DataVector& alpha,
                                     const DataVector& point) {
  typedef std::vector<std::pair<size_t, double> > IndexValVector;

  IndexValVector vec;
  LinearModifiedBasis<unsigned int, unsigned int> base;
  GetAffectedBasisFunctions <
  LinearModifiedBasis<unsigned int, unsigned int> > ga(storage);

  /* Scale point to bounding box */

  // Initialize a copy of point
  DataVector point_bb = DataVector(point.getSize());
  point_bb.copyFrom(point);

  // Get bounding box
  BoundingBox* bb = storage.getBoundingBox();
  size_t dim = bb->getDimension();

  if (bb != nullptr) {
    for (size_t d = 0; d < dim; ++d) {
      BoundingBox1D dimbb = bb->getBoundary(d);

      if (dimbb.leftBoundary == 0.0 && dimbb.rightBoundary == 1.0) {
        continue;
      }

      if (!(dimbb.leftBoundary <= point[d] &&
            point[d] <= dimbb.rightBoundary)) {
        return 0.0;
      }

      point_bb[d] = (point[d] - dimbb.leftBoundary) / (dimbb.rightBoundary -
                    dimbb.leftBoundary);
    }
  }

  ga(base, point_bb, vec);

  double result = 0.0;

  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}

}  // namespace base
}  // namespace sgpp
