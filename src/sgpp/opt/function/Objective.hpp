#ifndef SGPP_OPT_FUNCTION_OBJECTIVE_HPP
#define SGPP_OPT_FUNCTION_OBJECTIVE_HPP

#include <vector>
#include <cstddef>

namespace sg
{
namespace opt
{
namespace function
{

class Objective
{
public:
    Objective(size_t d) : d(d) {}
    virtual ~Objective() {}
    
    virtual double eval(const std::vector<double> &x) = 0;
    
    size_t getDimension() const { return d; }
    
protected:
    size_t d;
};

}
}
}

#endif
