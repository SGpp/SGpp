#ifndef SGPP_OPT_OPTIMIZATION_OPTIMIZERNELDERMEAD_HPP
#define SGPP_OPT_OPTIMIZATION_OPTIMIZERNELDERMEAD_HPP

#include "opt/optimization/Optimizer.hpp"

namespace sg
{
namespace opt
{
namespace optimization
{

class OptimizerNelderMead : public Optimizer
{
public:
    static const double DEFAULT_ALPHA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_RHO;
    static const double DEFAULT_SIGMA;
    static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 200;
    static const double STARTING_SIMPLEX_EDGE_LENGTH;
    
    OptimizerNelderMead(function::ObjectiveFunction &f);
    
    OptimizerNelderMead(function::ObjectiveFunction &f,
                        size_t max_it_count, double alpha, double gamma, double rho, double sigma,
                        size_t max_fcn_eval_count);
    
    void optimize(std::vector<double> &xopt);
    
    double getAlpha() const;
    void setAlpha(double alpha);
    
    double getGamma() const;
    void setGamma(double gamma);
    
    double getRho() const;
    void setRho(double rho);
    
    double getSigma() const;
    void setSigma(double sigma);
    
    size_t getMaxFcnEvalCount() const;
    void setMaxFcnEvalCount(size_t max_fcn_eval_count);
    
protected:
    double alpha;
    double gamma;
    double rho;
    double sigma;
    size_t max_fcn_eval_count;
};

}
}
}

#endif
