#ifndef SGPP_OPT_FUNCTION_OBJECTIVEGRADIENT_HPP
#define SGPP_OPT_FUNCTION_OBJECTIVEGRADIENT_HPP

#include "base/datatypes/DataVector.hpp"

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{   
namespace function
{

class ObjectiveGradient
{
public:
    ObjectiveGradient(size_t d) : d(d) {}
    virtual ~ObjectiveGradient() {}
    
    virtual double evalGradient(const std::vector<double> &x,
                                base::DataVector &gradient) = 0;
    
    size_t getDimensionCount() const { return d; }
    
protected:
    size_t d;
};

}
}
}

#endif
