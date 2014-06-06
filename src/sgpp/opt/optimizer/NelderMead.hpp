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
    static const double DEFAULT_BETA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_DELTA;
    static const size_t DEFAULT_MAX_IT_COUNT = 1000;
    static const double STARTING_SIMPLEX_EDGE_LENGTH;
    
    NelderMead(function::Objective &f,
               size_t max_it_count = DEFAULT_MAX_IT_COUNT,
               double alpha = DEFAULT_ALPHA,
               double beta = DEFAULT_BETA,
               double gamma = DEFAULT_GAMMA,
               double delta = DEFAULT_DELTA);
    
    double optimize(std::vector<double> &xopt);
    
    std::unique_ptr<Optimizer> clone();
    
    double getAlpha() const;
    void setAlpha(double alpha);
    
    double getBeta() const;
    void setBeta(double beta);
    
    double getGamma() const;
    void setGamma(double gamma);
    
    double getDelta() const;
    void setDelta(double delta);
    
    //size_t getMaxFcnEvalCount() const;
    //void setMaxFcnEvalCount(size_t max_fcn_eval_count);
    
protected:
    double alpha;
    double beta;
    double gamma;
    double delta;
    //size_t max_fcn_eval_count;
};

}
}
}

#endif
