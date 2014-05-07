#ifndef SGPP_OPT_OPTIMIZER_GRADIENTMETHOD_HPP
#define SGPP_OPT_OPTIMIZER_GRADIENTMETHOD_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveFunctionGradient.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

class GradientMethod : public Optimizer
{
public:
    static const size_t DEFAULT_MAX_IT_COUNT = 2000;
    static const double DEFAULT_BETA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_TOLERANCE;
    
    GradientMethod(function::ObjectiveFunction &f,
                   function::ObjectiveFunctionGradient &f_gradient);
    
    GradientMethod(function::ObjectiveFunction &f,
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
