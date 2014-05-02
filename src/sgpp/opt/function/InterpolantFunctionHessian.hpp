#ifndef SGPP_OPT_FUNCTION_INTERPOLANTFUNCTIONHESSIAN_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTFUNCTIONHESSIAN_HPP

#include "opt/function/ObjectiveFunctionHessian.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/OperationEvalHessian.hpp"

#include <vector>

namespace sg
{
namespace opt
{
namespace function
{

class InterpolantFunctionHessian : public ObjectiveFunctionHessian
{
public:
    InterpolantFunctionHessian(size_t d, base::OperationEvalHessian *op_eval_hessian,
                               base::DataVector &alpha) :
            ObjectiveFunctionHessian(d), op_eval_hessian(op_eval_hessian), alpha(alpha) {}
    
    inline double evalHessian(const std::vector<double> &x,
                              base::DataVector &gradient, base::DataMatrix &hessian)
    {
        std::vector<double> y = x;
        return op_eval_hessian->evalHessian(alpha, y, gradient, hessian);
    }
    
protected:
    base::OperationEvalHessian *op_eval_hessian;
    base::DataVector &alpha;
};

}
}
}

#endif
