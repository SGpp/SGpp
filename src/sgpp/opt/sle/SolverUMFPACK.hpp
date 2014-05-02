#ifndef SGPP_OPT_SLE_SOLVERUMFPACK_HPP
#define SGPP_OPT_SLE_SOLVERUMFPACK_HPP

#include "opt/sle/Solver.hpp"

#include <cstddef>
#include <vector>
#include <cstdint>

namespace sg
{
namespace opt
{
namespace sle
{

class SolverUMFPACK : public Solver
{
public:
    bool solve(System &system, std::vector<double> &x) const;
    bool solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
               const std::vector<double> &Tx, const std::vector<double> &b,
               std::vector<double> &x, bool initial_output = true) const;
};

}
}
}

#endif
