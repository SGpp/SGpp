#ifndef SGPP_OPT_FUNCTION_OBJECTIVEHESSIAN_HPP
#define SGPP_OPT_FUNCTION_OBJECTIVEHESSIAN_HPP

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

class ObjectiveHessian
{
public:
    ObjectiveHessian(size_t d) : d(d) {}
    virtual ~ObjectiveHessian() {}
    
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
