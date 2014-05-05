#ifndef SGPP_OPT_SLE_SYSTEM_PARALLELIZABLE_HPP
#define SGPP_OPT_SLE_SYSTEM_PARALLELIZABLE_HPP

#include "opt/sle/system/System.hpp"

namespace sg
{
namespace opt
{
namespace sle
{
namespace system
{

class Parallelizable : public System
{
public:
    Parallelizable(size_t n) : System(n) {}
    Parallelizable(const std::vector<double> &b) : System(b) {}
    virtual ~Parallelizable() {}
    
    virtual Parallelizable *clone() = 0;
    
    bool isParallelizable() const { return true; }
};

}
}
}
}

#endif
