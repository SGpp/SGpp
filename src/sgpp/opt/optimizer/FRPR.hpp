#ifndef SGPP_OPT_OPTIMIZER_FRPR_HPP
#define SGPP_OPT_OPTIMIZER_FRPR_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveFunctionGradient.hpp"

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
    
    FRPR(function::ObjectiveFunction &f, function::ObjectiveFunctionGradient &f_gradient);
    
    FRPR(function::ObjectiveFunction &f,
         function::ObjectiveFunctionGradient &f_gradient,
         size_t max_it_count, double tolerance);
    
    void optimize(std::vector<double> &xopt);
    
    function::ObjectiveFunctionGradient &getObjectiveFunctionGradient() const;
    
    double getTolerance() const;
    void setTolerance(double tolerance);
    
protected:
    function::ObjectiveFunctionGradient &f_gradient;
    double tol;
};

}
}
}

#endif
