#include "opt/optimization/OptimizerFRPR.hpp"
#include "opt/optimization/numrecipes/ObjectiveFunction.hpp"
#include "opt/optimization/numrecipes/Frprmn.hpp"
#include "opt/tools/Printer.hpp"

namespace sg
{
namespace opt
{
namespace optimization
{

const double OptimizerFRPR::DEFAULT_TOLERANCE = 1e-20;

OptimizerFRPR::OptimizerFRPR(
        function::ObjectiveFunction &f,
        function::ObjectiveFunctionGradient &f_gradient) :
    OptimizerFRPR(f, f_gradient, DEFAULT_MAX_IT_COUNT, DEFAULT_TOLERANCE)
{
}

OptimizerFRPR::OptimizerFRPR(function::ObjectiveFunction &f,
                             function::ObjectiveFunctionGradient &f_gradient,
                             size_t max_it_count, double tolerance) :
    Optimizer(f, max_it_count),
    f_gradient(f_gradient),
    tol(tolerance)
{
}

void OptimizerFRPR::optimize(std::vector<double> &xopt)
{
    tools::printer.printStatusBegin("Optimizing (FRPR)...");
    
    numrecipes::ObjectiveFunction num_rec_obj_fcn(f, f_gradient);
    numrecipes::Frprmn<numrecipes::ObjectiveFunction> frprmn(num_rec_obj_fcn, tol, N);
    xopt = frprmn.minimize(x0);
    
    tools::printer.printStatusEnd();
}

function::ObjectiveFunctionGradient &OptimizerFRPR::getObjectiveFunctionGradient() const
{
    return f_gradient;
}

double OptimizerFRPR::getTolerance() const
{
    return tol;
}

void OptimizerFRPR::setTolerance(double tolerance)
{
    tol = tolerance;
}

}
}
}
