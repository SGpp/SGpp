#ifndef SGPP_OPT_FUNCTION_OBJECTIVEFUNCTIONHESSIAN_HPP
#define SGPP_OPT_FUNCTION_OBJECTIVEFUNCTIONHESSIAN_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace function
{

class ObjectiveFunctionHessian
{
public:
    ObjectiveFunctionHessian(size_t d) : d(d) {}
    virtual ~ObjectiveFunctionHessian() {}
    
    virtual double evalHessian(const std::vector<double> &x,
                               base::DataVector &gradient, base::DataMatrix &hessian) = 0;
    
    size_t getDimensionCount() const { return d; }
    
protected:
    size_t d;
};

}
}
}

#endif
