#ifndef SGPP_OPT_FUNCTION_INTERPOLANTHESSIAN_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTHESSIAN_HPP

#include "opt/function/ObjectiveHessian.hpp"
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

class InterpolantHessian : public ObjectiveHessian
{
public:
    InterpolantHessian(size_t d, base::Grid &grid, base::DataVector &alpha) :
            ObjectiveHessian(d), grid(grid),
            op_eval_hessian(op_factory::createOperationEvalHessian(grid)),
            alpha(alpha) {}
    
    inline double evalHessian(const std::vector<double> &x,
                              base::DataVector &gradient, base::DataMatrix &hessian)
    {
        std::vector<double> y = x;
        return op_eval_hessian->evalHessian(alpha, y, gradient, hessian);
    }
    
    virtual std::unique_ptr<ObjectiveHessian> clone()
    {
        return std::unique_ptr<ObjectiveHessian>(new InterpolantHessian(d, grid, alpha));
    }
    
protected:
    base::Grid &grid;
    std::unique_ptr<base::OperationEvalHessian> op_eval_hessian;
    base::DataVector &alpha;
};

}
}
}

#endif
