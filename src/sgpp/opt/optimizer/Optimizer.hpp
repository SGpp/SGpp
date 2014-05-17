#ifndef SGPP_OPT_OPTIMIZER_OPTIMIZER_HPP
#define SGPP_OPT_OPTIMIZER_OPTIMIZER_HPP

#include "opt/function/Objective.hpp"

#include <vector>
#include <cstddef>

#include <iostream>

namespace sg
{
namespace opt
{
namespace optimizer
{

class Optimizer
{
public:
    static const size_t DEFAULT_MAX_IT_COUNT = 200;
    
    Optimizer(function::Objective &f, size_t max_it_count = DEFAULT_MAX_IT_COUNT) :
            f(f.clone()), N(max_it_count), x0(std::vector<double>(f.getDimension(), 0.5)) {}
    
    Optimizer(Optimizer &&other) : f(std::move(other.f)), N(other.N), x0(other.x0) {}
    
    virtual ~Optimizer() {}
    
    virtual double optimize(std::vector<double> &xopt) = 0;
    virtual std::unique_ptr<Optimizer> clone() = 0;
    
    const std::unique_ptr<function::Objective> &getObjectiveFunction() const { return f; }
    
    size_t getMaxItCount() const { return N; }
    void setMaxItCount(size_t max_it_count) { N = max_it_count; }
    
    const std::vector<double> &getStartingPoint() const { return x0; }
    void setStartingPoint(const std::vector<double> &x0) { this->x0 = x0; }
    
protected:
    std::unique_ptr<function::Objective> f;
    size_t N;
    std::vector<double> x0;
};

}
}
}

#endif
