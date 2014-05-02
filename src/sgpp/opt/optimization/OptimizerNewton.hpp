#ifndef SGPP_OPT_OPTIMIZATION_OPTIMIZERNEWTON_HPP
#define SGPP_OPT_OPTIMIZATION_OPTIMIZERNEWTON_HPP

#include "opt/optimization/Optimizer.hpp"
#include "opt/function/ObjectiveFunctionHessian.hpp"
#include "opt/sle/Solver.hpp"
#include "opt/sle/SolverBiCGStab.hpp"

#include <cstddef>

namespace sg
{
namespace opt
{
namespace optimization
{

class OptimizerNewton : public Optimizer
{
public:
    static const double DEFAULT_ALPHA1;
    static const double DEFAULT_ALPHA2;
    static const double DEFAULT_BETA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_P;
    static const double DEFAULT_TOLERANCE;
    
    OptimizerNewton(function::ObjectiveFunction &f,
                    function::ObjectiveFunctionHessian &f_hessian);
    
    OptimizerNewton(function::ObjectiveFunction &f,
                    function::ObjectiveFunctionHessian &f_hessian,
                    size_t max_it_count, double alpha1, double alpha2, double beta, double gamma,
                    double p, double tolerance, const sle::Solver &sle_solver);
    
    void optimize(std::vector<double> &xopt);
    
    function::ObjectiveFunctionHessian &getObjectiveFunctionHessian() const;
    
    double getAlpha1() const;
    void setAlpha1(double alpha1);
    
    double getAlpha2() const;
    void setAlpha2(double alpha1);
    
    double getBeta() const;
    void setBeta(double beta);
    
    double getGamma() const;
    void setGamma(double gamma);
    
    double getP() const;
    void setP(double p);
    
    double getTolerance() const;
    void setTolerance(double tolerance);
    
protected:
    function::ObjectiveFunctionHessian &f_hessian;
    double alpha1;
    double alpha2;
    double beta;
    double gamma;
    double p;
    double tol;
    const sle::SolverBiCGStab default_sle_solver;
    const sle::Solver &sle_solver;
};

}
}
}

#endif
