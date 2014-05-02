#ifndef SGPP_OPT_FUNCTION_OBJECTIVEFUNCTION_HPP
#define SGPP_OPT_FUNCTION_OBJECTIVEFUNCTION_HPP

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace function
{

class ObjectiveFunction
{
public:
    ObjectiveFunction(size_t d) : d(d) {}
    virtual ~ObjectiveFunction() {}
    
    virtual double eval(const std::vector<double> &x) = 0;
    
    size_t getDimension() const { return d; }
    
protected:
    size_t d;
};

}
}
}

#endif
