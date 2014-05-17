#ifndef SGPP_OPT_OPTIMIZER_NEWTON_HPP
#define SGPP_OPT_OPTIMIZER_NEWTON_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/function/ObjectiveHessian.hpp"
#include "opt/sle/solver/Solver.hpp"
#include "opt/sle/solver/BiCGStab.hpp"

#include <cstddef>

namespace sg
{
namespace opt
{
namespace optimizer
{

class Newton : public Optimizer
{
public:
    static const double DEFAULT_ALPHA1;
    static const double DEFAULT_ALPHA2;
    static const double DEFAULT_BETA;
    static const double DEFAULT_GAMMA;
    static const double DEFAULT_P;
    static const double DEFAULT_TOLERANCE;
    
    Newton(function::Objective &f, function::ObjectiveHessian &f_hessian);
    
    Newton(function::Objective &f, function::ObjectiveHessian &f_hessian,
           size_t max_it_count, double alpha1, double alpha2, double beta, double gamma,
           double p, double tolerance, const sle::solver::Solver &sle_solver);
    
    double optimize(std::vector<double> &xopt);
    
    std::unique_ptr<Optimizer> clone();
    
    const std::unique_ptr<function::ObjectiveHessian> &getObjectiveHessian() const;
    
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
    std::unique_ptr<function::ObjectiveHessian> f_hessian;
    double alpha1;
    double alpha2;
    double beta;
    double gamma;
    double p;
    double tol;
    const sle::solver::BiCGStab default_sle_solver;
    const sle::solver::Solver &sle_solver;
    
    void initialize(double alpha1, double alpha2, double beta, double gamma,
                    double p, double tolerance);
};

}
}
}

#endif
