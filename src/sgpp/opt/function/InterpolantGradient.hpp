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
    InterpolantGradient(size_t d, base::Grid &grid, base::DataVector &alpha) :
            ObjectiveGradient(d), grid(grid),
            op_eval_gradient(op_factory::createOperationEvalGradient(grid)),
            alpha(alpha) {}
    
    inline double evalGradient(const std::vector<double> &x, base::DataVector &gradient)
    {
        std::vector<double> y = x;
        return op_eval_gradient->evalGradient(alpha, y, gradient);
    }
    
    virtual std::unique_ptr<ObjectiveGradient> clone()
    {
        return std::unique_ptr<ObjectiveGradient>(new InterpolantGradient(d, grid, alpha));
    }
    
protected:
    base::Grid &grid;
    std::unique_ptr<base::OperationEvalGradient> op_eval_gradient;
    base::DataVector &alpha;
};

}
}
}

#endif
