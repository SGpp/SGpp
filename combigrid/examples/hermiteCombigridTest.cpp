#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/PsiHermiteInterpolationEvaluator.hpp>
#include <sgpp_combigrid.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

double f(sgpp::base::DataVector const &x) {
  
  return x[0];
  
}

// We have to wrap f in a sgpp::combigrid::MultiFunction object.
static const sgpp::combigrid::MultiFunction func(f);

int main() {

/*
  std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy> hierarchy =
      sgpp::combigrid::CombiHierarchies::expUniform();

  std::shared_ptr<sgpp::combigrid::PsiHermiteInterpolationEvaluator> evaluator =
      std::make_shared<sgpp::combigrid::PsiHermiteInterpolationEvaluator>();

  std::vector<double> gridpoints = hierarchy->getPoints(3, true);

  std::vector<double> alpha(gridpoints.size(), 0.0);
*/
  size_t d = 1;
std::shared_ptr<sgpp::combigrid::CombigridOperation> operation =
      sgpp::combigrid::CombigridOperation::createExpUniformBoundaryPsiLinearInterpolation(d,0, func);


sgpp::base::DataVector evaluationPoint(d);

  evaluationPoint[0] = 0.4;

 



    double result = operation->evaluate(1, evaluationPoint);

 
  std::cout << "Interpolation result: " << result << ", function value: " << func(evaluationPoint)
            << "\n";
/*
  for (size_t i = 0; i < gridpoints.size(); i++) {
    std::cout << gridpoints[i] << std::endl;
    alpha[i] = f(gridpoints[i]);
  }

  sgpp::combigrid::FloatScalarVector evalpoint = sgpp::combigrid::FloatScalarVector(0.5);

  evaluator->setGridPoints(gridpoints);
  evaluator->setParameter(evalpoint);

  sgpp::combigrid::FloatScalarVector y = evaluator->eval(alpha);

  std::cout << y.value();

  */
}