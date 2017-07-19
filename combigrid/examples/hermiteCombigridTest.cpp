#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/PsiHermiteInterpolationEvaluator.hpp>
#include <sgpp_combigrid.hpp>

double f(double x) { return 1; }

int main() {
  int numDimensions = 1;

  std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy> hierarchy =
      sgpp::combigrid::CombiHierarchies::expUniform();

  std::shared_ptr<sgpp::combigrid::PsiHermiteInterpolationEvaluator> evaluator =
      std::make_shared<sgpp::combigrid::PsiHermiteInterpolationEvaluator>();

  std::vector<double> gridpoints = hierarchy->getPoints(3, true);

  std::vector<double> alpha(gridpoints.size(), 0.0);

  for (size_t i = 0; i < gridpoints.size(); i++) {
    std::cout << gridpoints[i] << std::endl;
    alpha[i] = f(gridpoints[i]);
  }

  sgpp::combigrid::FloatScalarVector evalpoint = sgpp::combigrid::FloatScalarVector(0.5);

  evaluator->setGridPoints(gridpoints);
  evaluator->setParameter(evalpoint);

  sgpp::combigrid::FloatScalarVector y = evaluator->eval(alpha);

  std::cout << y.value();
}