#ifndef SGPP_OPT_OPTIMIZATION_OPTIMIZERFRPR_HPP
#define SGPP_OPT_OPTIMIZATION_OPTIMIZERFRPR_HPP

#include "opt/optimization/Optimizer.hpp"
#include "opt/function/ObjectiveFunctionGradient.hpp"

namespace sg
{
namespace opt
{
namespace optimization
{

class OptimizerFRPR : public Optimizer
{
public:
    static const double DEFAULT_TOLERANCE;
    
    OptimizerFRPR(function::ObjectiveFunction &f,
                  function::ObjectiveFunctionGradient &f_gradient);
    
    OptimizerFRPR(function::ObjectiveFunction &f,
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
