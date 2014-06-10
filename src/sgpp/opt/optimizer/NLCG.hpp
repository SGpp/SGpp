#ifndef SGPP_OPT_OPTIMIZER_NLCG_HPP
#define SGPP_OPT_OPTIMIZER_NLCG_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveGradient.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

class NLCG : public Optimizer
{
public:
    static const double DEFAULT_BETA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_TOLERANCE;
    static const double DEFAULT_EPSILON;
    static const double DEFAULT_RESTART_THRESHOLD;
    
    NLCG(function::Objective &f,
         function::ObjectiveGradient &f_gradient,
         size_t max_it_count = DEFAULT_MAX_IT_COUNT,
         double beta = DEFAULT_BETA,
         double gamma = DEFAULT_GAMMA,
         double tolerance = DEFAULT_TOLERANCE,
         double epsilon = DEFAULT_EPSILON,
         double restart_threshold = DEFAULT_RESTART_THRESHOLD);
    
    double optimize(std::vector<double> &xopt);
    
    std::unique_ptr<Optimizer> clone();
    
    const std::unique_ptr<function::ObjectiveGradient> &getObjectiveGradient() const;
    
    double getBeta() const;
    void setBeta(double beta);
    
    double getGamma() const;
    void setGamma(double gamma);
    
    double getTolerance() const;
    void setTolerance(double tolerance);
    
    double getEpsilon() const;
    void setEpsilon(double epsilon);
    
    double getRestartThreshold() const;
    void setRestartThreshold(double restart_threshold);
    
protected:
    std::unique_ptr<function::ObjectiveGradient> f_gradient;
    double beta;
    double gamma;
    double tol;
    double eps;
    double alpha;
};

}
}
}

#endif
