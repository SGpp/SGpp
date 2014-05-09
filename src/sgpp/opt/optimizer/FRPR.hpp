#ifndef SGPP_OPT_OPTIMIZER_FRPR_HPP
#define SGPP_OPT_OPTIMIZER_FRPR_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveGradient.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

class FRPR : public Optimizer
{
public:
    static const double DEFAULT_TOLERANCE;
    
    FRPR(function::Objective &f, function::ObjectiveGradient &f_gradient);
    
    FRPR(function::Objective &f, function::ObjectiveGradient &f_gradient,
         size_t max_it_count, double tolerance);
    
    void optimize(std::vector<double> &xopt);
    
    function::ObjectiveGradient &getObjectiveGradient() const;
    
    double getTolerance() const;
    void setTolerance(double tolerance);
    
protected:
    function::ObjectiveGradient &f_gradient;
    double tol;
};

}
}
}

#endif
