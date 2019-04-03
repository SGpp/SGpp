// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_pce_cpp PCE with Combigrids
 *
 * This simple example shows how to create a Polynomial Chaos Expansion from an
 * adaptively refined combigrid.
 */

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <cmath>

#include <iostream>
#include <vector>

int main() {
  try {
    /**
     * First we have to define a model to approximate.
     */
    sgpp::combigrid::Ishigami ishigamiModel;
    sgpp::combigrid::MultiFunction func(ishigamiModel.eval);

    /**
     *  Then we can create a refined combigrid
     */

    // create polynomial basis
    sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
    config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
    auto basisFunction = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

    for (size_t q = 6; q < 7; ++q) {
      // create sprarse grid interpolation operation
      auto tensor_op =
          sgpp::combigrid::CombigridTensorOperation::createExpClenshawCurtisPolynomialInterpolation(
              basisFunction, ishigamiModel.numDims, func);
      sgpp::combigrid::Stopwatch stopwatch;
      stopwatch.start();
      // start with regular level q and add some level adaptively
      tensor_op->getLevelManager()->addRegularLevels(q);
      tensor_op->getLevelManager()->addLevelsAdaptiveByNumLevels(10);
      tensor_op->getLevelManager()->addLevelsAdaptiveByNumLevels(10);
      stopwatch.log();
      stopwatch.start();

      /**
       * and construct a PCE representation to easily calculate statistical features of our model.
       */

      // create polynomial chaos surrogate from sparse grid
      sgpp::combigrid::CombigridSurrogateModelConfiguration config;
      config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
      config.loadFromCombigridOperation(tensor_op, false);
      config.basisFunction = basisFunction;
      auto pce = sgpp::combigrid::createCombigridSurrogateModel(config);

      stopwatch.log();
      stopwatch.start();

      // compute mean, variance and sobol indices
      double mean = pce->mean();
      double variance = pce->variance();
      sgpp::base::DataVector sobol_indices;
      sgpp::base::DataVector total_sobol_indices;
      pce->getComponentSobolIndices(sobol_indices);
      pce->getTotalSobolIndices(total_sobol_indices);

      // print results
      std::cout << "Time: "
                << stopwatch.elapsedSeconds() / static_cast<double>(tensor_op->numGridPoints())
                << std::endl;
      std::cout << "---------------------------------------------------------" << std::endl;
      std::cout << "#gp = " << tensor_op->getLevelManager()->numGridPoints() << std::endl;
      std::cout << "E(u) = " << mean << std::endl;
      std::cout << "Var(u) = " << variance << std::endl;
      std::cout << "Sobol indices = " << sobol_indices.toString() << std::endl;
      std::cout << "Sum Sobol indices = " << sobol_indices.sum() << std::endl;
      std::cout << "Total Sobol indices = " << total_sobol_indices.toString() << std::endl;
      std::cout << "---------------------------------------------------------" << std::endl;
    }
  }
  catch (sgpp::base::generation_exception& exc)  {
    std::cout << "Exception: " << exc.what() << std::endl;
    std::cout << "Skipping example..." << std::endl;
  }
}

/**
* Output:
* @verbatim
* Time: 3.74635s.
* Time: 3.73112s.
* Time: 4.96569e-05
* ---------------------------------------------------------
* #gp = 1825
* E(u) = 3.5
* Var(u) = 13.8446
* Sobol indices = [3.13905191147811180041e-01, 4.42411144790040733454e-01,
* 9.56029390935037928152e-31, 5.58403133152096916677e-32, 2.43683664062148142015e-01,
* 6.16110611418130297137e-32, 8.27081513075557474542e-32]
* Sum Sobol indices = 1
* Total Sobol indices = [5.57588855209959377568e-01, 4.42411144790040733454e-01,
* 2.43683664062148142015e-01]
* @endverbatim
*/
