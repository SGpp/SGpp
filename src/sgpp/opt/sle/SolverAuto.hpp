#ifndef SGPP_OPT_SLE_SOLVERAUTO_HPP
#define SGPP_OPT_SLE_SOLVERAUTO_HPP

#include "opt/sle/Solver.hpp"

#include <vector>
#include <cstdint>

namespace sg
{
namespace opt
{
namespace sle
{

class SolverAuto : public Solver
{
public:
    static const double MAX_NNZ_RATIO_FOR_SPARSE;
    static const double ESTIMATE_NNZ_ROWS_SAMPLE_SIZE;
    
    bool solve(System &system, std::vector<double> &x) const;
    /*bool solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
               const std::vector<double> &Tx, const std::vector<double> &b,
               std::vector<double> &x) const;*/
};

}
}
}

#endif
