#ifndef SGPP_OPT_OPTIMIZATION_OPTIMIZERGRADIENTMETHOD_HPP
#define SGPP_OPT_OPTIMIZATION_OPTIMIZERGRADIENTMETHOD_HPP

#include "opt/optimization/Optimizer.hpp"
#include "opt/function/ObjectiveFunctionGradient.hpp"

namespace sg
{
namespace opt
{
namespace optimization
{

class OptimizerGradientMethod : public Optimizer
{
public:
    static const size_t DEFAULT_MAX_IT_COUNT = 2000;
    static const double DEFAULT_BETA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_TOLERANCE;
    
    OptimizerGradientMethod(function::ObjectiveFunction &f,
                            function::ObjectiveFunctionGradient &f_gradient);
    
    OptimizerGradientMethod(function::ObjectiveFunction &f,
                            function::ObjectiveFunctionGradient &f_gradient,
                            size_t max_it_count, double beta, double gamma, double tolerance);
    
    void optimize(std::vector<double> &xopt);
    
    function::ObjectiveFunctionGradient &getObjectiveFunctionGradient() const;
    
    double getBeta() const;
    void setBeta(double beta);
    
    double getGamma() const;
    void setGamma(double gamma);
    
    double getTolerance() const;
    void setTolerance(double tolerance);
    
protected:
    function::ObjectiveFunctionGradient &f_gradient;
    double beta;
    double gamma;
    double tol;
};

}
}
}

#endif
