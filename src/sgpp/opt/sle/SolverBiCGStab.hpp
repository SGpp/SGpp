#ifndef SGPP_OPT_SLE_SOLVERBICGSTAB_HPP
#define SGPP_OPT_SLE_SOLVERBICGSTAB_HPP

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

class SolverBiCGStab : public Solver
{
public:
    static const size_t DEFAULT_MAX_IT_COUNT = 1000;
    static const double DEFAULT_TOLERANCE;
    
    SolverBiCGStab();
    SolverBiCGStab(size_t max_it_count, double tolerance, const std::vector<double> &x0);
    
    bool solve(System &system, std::vector<double> &x) const;
    /*bool solve(const std::vector<uint32_t> &Ti, const std::vector<uint32_t> &Tj,
               const std::vector<double> &Tx, const std::vector<double> &b,
               std::vector<double> &x) const;*/
    
    size_t getMaxItCount() const;
    void setMaxItCount(size_t max_it_count);
    
    double getTolerance() const;
    void setTolerance(double tolerance);
    
    const std::vector<double> &getX0() const;
    void setX0(const std::vector<double> &x0);
    
protected:
    size_t N;
    double tol;
    std::vector<double> x0;
};

}
}
}

#endif
