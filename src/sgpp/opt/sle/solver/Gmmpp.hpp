#ifndef SGPP_OPT_SLE_SOLVER_GMMPP_HPP
#define SGPP_OPT_SLE_SOLVER_GMMPP_HPP

#include "opt/sle/solver/Solver.hpp"

#include <vector>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

class Gmmpp : public Solver
{
public:
    bool solve(system::System &system, std::vector<double> &x) const;
};

}
}
}
}

#endif
