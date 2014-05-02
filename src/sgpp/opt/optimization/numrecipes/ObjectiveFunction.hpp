#ifndef SGPP_OPT_OPTIMIZATION_NUMRECIPES_OBJECTIVEFUNCTIONGRADIENT_HPP
#define SGPP_OPT_OPTIMIZATION_NUMRECIPES_OBJECTIVEFUNCTIONGRADIENT_HPP

#include "base/datatypes/DataVector.hpp"
#include "opt/function/ObjectiveFunction.hpp"
#include "opt/function/ObjectiveFunctionGradient.hpp"

#include <vector>
#include <cstddef>
#include <cmath>

namespace sg
{
namespace opt
{
namespace optimization
{
namespace numrecipes
{

struct ObjectiveFunction
{
    function::ObjectiveFunction &f;
    function::ObjectiveFunctionGradient &f_gradient;
    
    ObjectiveFunction(function::ObjectiveFunction &f,
                      function::ObjectiveFunctionGradient &f_gradient) :
        f(f),
        f_gradient(f_gradient)
    {
    }
    
    double operator()(const std::vector<double> &x)
    {
        for (size_t i = 0; i < x.size(); i++)
        {
            if ((x[i] < 0.0) || (x[i] > 1.0))
            {
                return INFINITY;
            }
        }
        
        return f.eval(x);
    }
    
    void df(const std::vector<double> &x, std::vector<double> &grad_f)
    {
        base::DataVector grad_f_datavec(x.size());
        f_gradient.evalGradient(x, grad_f_datavec);
        grad_f = std::vector<double>(grad_f_datavec.getPointer(),
                                     grad_f_datavec.getPointer() + grad_f_datavec.getSize());
    }
};

}
}
}
}

#endif
