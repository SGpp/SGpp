#ifndef SGPP_OPT_SLE_SOLVER_HPP
#define SGPP_OPT_SLE_SOLVER_HPP

#include "opt/sle/System.hpp"

#include <vector>
#include <cstdint>

namespace sg
{
namespace opt
{
namespace sle
{

class Solver
{
public:
    Solver() {}
    virtual ~Solver() {}
    
    virtual bool solve(System &system, std::vector<double> &x) const = 0;
    /*virtual bool solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
                       const std::vector<double> &Tx, const std::vector<double> &b,
                       std::vector<double> &x) const = 0;*/
};

}
}
}

#endif
