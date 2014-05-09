#ifndef SGPP_OPT_FUNCTION_INTERPOLANTGRADIENT_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTGRADIENT_HPP

#include "opt/function/ObjectiveGradient.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationEvalGradient.hpp"

#include <vector>

namespace sg
{
namespace opt
{
namespace function
{

class InterpolantGradient : public ObjectiveGradient
{
public:
    InterpolantGradient(size_t d, base::OperationEvalGradient *op_eval_gradient,
                                base::DataVector &alpha) :
            ObjectiveGradient(d), op_eval_gradient(op_eval_gradient), alpha(alpha) {}
    
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
