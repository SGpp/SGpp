#ifndef SGPP_OPT_FUNCTION_INTERPOLANTFUNCTIONGRADIENT_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTFUNCTIONGRADIENT_HPP

#include "opt/function/ObjectiveFunctionGradient.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationEvalGradient.hpp"

#include <vector>

namespace sg
{
namespace opt
{
namespace function
{

class InterpolantFunctionGradient : public ObjectiveFunctionGradient
{
public:
    InterpolantFunctionGradient(size_t d, base::OperationEvalGradient *op_eval_gradient,
                                base::DataVector &alpha) :
            ObjectiveFunctionGradient(d), op_eval_gradient(op_eval_gradient), alpha(alpha) {}
    
    inline double evalGradient(const std::vector<double> &x, base::DataVector &gradient)
    {
        std::vector<double> y = x;
        return op_eval_gradient->evalGradient(alpha, y, gradient);
    }
    
protected:
    base::OperationEvalGradient *op_eval_gradient;
    base::DataVector &alpha;
};

}
}
}

#endif
