#ifndef SGPP_OPT_SLE_SOLVER_AUTO_HPP
#define SGPP_OPT_SLE_SOLVER_AUTO_HPP

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

class Auto : public Solver
{
public:
    static const size_t MAX_DIM_FOR_FULL = 30000;
    static const double MAX_NNZ_RATIO_FOR_SPARSE;
    static const double ESTIMATE_NNZ_ROWS_SAMPLE_SIZE;
    static const double MAX_NNZ_RATIO_FOR_GMMPP;
    
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
