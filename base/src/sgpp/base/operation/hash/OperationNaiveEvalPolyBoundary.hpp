/*
 * OperationNaiveEvalPolyBoundary.hpp
 *
 *  Created on: Jul 1, 2015
 *      Author: franzefn
 */

#ifndef OPERATIONNAIVEEVALPOLYBOUNDARY_HPP_
#define OPERATIONNAIVEEVALPOLYBOUNDARY_HPP_

#include "OperationNaiveEval.hpp"

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace base {

class OperationNaiveEvalPolyBoundary: public base::OperationNaiveEval {
public:

    /**
     * Constructor.
     *
     * @param storage   storage of the sparse grid
     * @param degree    polynomial degree
     */
    OperationNaiveEvalPolyBoundary(base::GridStorage* storage, size_t degree) :
            storage(storage), base(degree) {
    }

    virtual ~OperationNaiveEvalPolyBoundary() {
    }

    /**
     * @param alpha     coefficient vector
     * @param point     evaluation point
     * @return          value of linear combination
     */
    virtual float_t eval(base::DataVector& alpha, base::DataVector& point);

protected:
    /// storage of the sparse grid
    base::GridStorage* storage;
    /// 1D B-spline basis
    SPolyBoundaryBase base;
};

} /* namespace base */
} /* namespace SGPP */

#endif /* OPERATIONNAIVEEVALPOLYBOUNDARY_HPP_ */
