#ifndef SGPP_OPT_SLE_SYSTEM_CLONEABLE_HPP
#define SGPP_OPT_SLE_SYSTEM_CLONEABLE_HPP

#include "opt/sle/system/System.hpp"

#include <memory>

namespace sg
{
namespace opt
{
namespace sle
{
namespace system
{

class Cloneable : public System
{
public:
    Cloneable(size_t n) : System(n) {}
    Cloneable(const std::vector<double> &b) : System(b) {}
    virtual ~Cloneable() {}
    
    virtual std::unique_ptr<Cloneable> clone() = 0;
    
    bool isCloneable() const { return true; }
};

}
}
}
}

#endif
