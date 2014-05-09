#include "opt/optimizer/FRPR.hpp"
#include "opt/optimizer/numrecipes/ObjectiveFunctionGradient.hpp"
#include "opt/optimizer/numrecipes/Frprmn.hpp"
#include "opt/tools/Printer.hpp"

namespace sg
{
namespace opt
{
namespace optimizer
{

const double FRPR::DEFAULT_TOLERANCE = 1e-20;

FRPR::FRPR(function::Objective &f, function::ObjectiveGradient &f_gradient) :
    FRPR(f, f_gradient, DEFAULT_MAX_IT_COUNT, DEFAULT_TOLERANCE)
{
}

FRPR::FRPR(function::Objective &f, function::ObjectiveGradient &f_gradient,
           size_t max_it_count, double tolerance) :
    Optimizer(f, max_it_count),
    f_gradient(f_gradient),
    tol(tolerance)
{
}

void FRPR::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (FRPR)...");
    
    numrecipes::ObjectiveFunctionGradient num_rec_obj_fcn(f, f_gradient);
    numrecipes::Frprmn<numrecipes::ObjectiveFunctionGradient> frprmn(num_rec_obj_fcn, tol, N);
    xopt = frprmn.minimize(x0);
    
    tools::printer.printStatusEnd();
}

function::ObjectiveGradient &FRPR::getObjectiveGradient() const
{
    return f_gradient;
}

double FRPR::getTolerance() const
{
    return tol;
}

void FRPR::setTolerance(double tolerance)
{
    tol = tolerance;
}

}
}
}
