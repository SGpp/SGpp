#ifndef SGPP_OPT_SLE_SOLVER_EIGEN_HPP
#define SGPP_OPT_SLE_SOLVER_EIGEN_HPP

#include "opt/sle/solver/Solver.hpp"

#include <vector>
#include <cstdint>

namespace sg
{
namespace opt
{
namespace sle
{
namespace solver
{

class Eigen : public Solver
{
public:
    bool solve(system::System &system, std::vector<double> &x) const;
    /*bool solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
               const std::vector<double> &Tx, const std::vector<double> &b,
               std::vector<double> &x) const;*/
};

}
}
}
}

#endif
