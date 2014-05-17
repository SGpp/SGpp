#ifndef SGPP_OPT_OPTIMIZER_NELDERMEAD_HPP
#define SGPP_OPT_OPTIMIZER_NELDERMEAD_HPP

#include "opt/optimizer/Optimizer.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

class NelderMead : public Optimizer
{
public:
    static const double DEFAULT_ALPHA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_RHO;
    static const double DEFAULT_SIGMA;
    static const size_t DEFAULT_MAX_IT_COUNT = 1000;
    static const double STARTING_SIMPLEX_EDGE_LENGTH;
    
    NelderMead(function::Objective &f,
               size_t max_it_count = DEFAULT_MAX_IT_COUNT,
               double alpha = DEFAULT_ALPHA,
               double gamma = DEFAULT_GAMMA,
               double rho = DEFAULT_RHO,
               double sigma = DEFAULT_SIGMA);
    
    double optimize(std::vector<double> &xopt);
    
    std::unique_ptr<Optimizer> clone();
    
    double getAlpha() const;
    void setAlpha(double alpha);
    
    double getGamma() const;
    void setGamma(double gamma);
    
    double getRho() const;
    void setRho(double rho);
    
    double getSigma() const;
    void setSigma(double sigma);
    
    //size_t getMaxFcnEvalCount() const;
    //void setMaxFcnEvalCount(size_t max_fcn_eval_count);
    
protected:
    double alpha;
    double gamma;
    double rho;
    double sigma;
    //size_t max_fcn_eval_count;
};

}
}
}

#endif
