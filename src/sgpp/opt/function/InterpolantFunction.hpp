#ifndef SGPP_OPT_FUNCTION_INTERPOLANTFUNCTION_HPP
#define SGPP_OPT_FUNCTION_INTERPOLANTFUNCTION_HPP

#include "opt/function/ObjectiveFunction.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationEval.hpp"

#include <vector>

namespace sg
{
namespace opt
{
namespace function
{

class InterpolantFunction : public ObjectiveFunction
{
public:
    InterpolantFunction(size_t d, base::OperationEval *op_eval, base::DataVector &alpha) :
            ObjectiveFunction(d), op_eval(op_eval), alpha(alpha) {}
    
    inline double eval(const std::vector<double> &x)
    {
        std::vector<double> y = x;
        return op_eval->eval(alpha, y);
    }
    
protected:
    base::OperationEval *op_eval;
    base::DataVector &alpha;
};

}
}
}

#endif
