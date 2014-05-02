#ifndef SGPP_OPT_SLE_SOLVEREIGEN_HPP
#define SGPP_OPT_SLE_SOLVEREIGEN_HPP

#include "opt/sle/Solver.hpp"

#include <vector>
#include <cstdint>

namespace sg
{
namespace opt
{
namespace sle
{

class SolverEigen : public Solver
{
public:
    bool solve(System &system, std::vector<double> &x) const;
    /*bool solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
               const std::vector<double> &Tx, const std::vector<double> &b,
               std::vector<double> &x) const;*/
};

}
}
}

#endif
