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
    //static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 200;
    static const size_t DEFAULT_MAX_IT_COUNT = 1000;
    static const double STARTING_SIMPLEX_EDGE_LENGTH;
    
    NelderMead(function::Objective &f);
    NelderMead(function::Objective &f, size_t max_it_count);
    /*NelderMead(function::Objective &f,
               size_t max_it_count, double alpha, double gamma, double rho, double sigma,
               size_t max_fcn_eval_count);*/
    NelderMead(function::Objective &f,
               size_t max_it_count, double alpha, double gamma, double rho, double sigma);
    
    void optimize(std::vector<double> &xopt);
    
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
