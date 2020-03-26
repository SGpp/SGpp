// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>

#include <string>

namespace sgpp {
namespace datadriven {
class DBMatPermutationFactory {
 public:
  DBMatPermutationFactory();
  explicit DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store);
  explicit DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store,
                                   const std::string& dbFilePath);

  /**
   * @brief Returns a offline object matching the specified configuration.
   * If the store contains a suitable base object, the permutation and blow-up approach is applied
   * on a copy of the base object, which is then returned.
   * If no suitable base object exists in the store, the factory builds a suitable object from
   * scratch and stores it. If a database file is given, the database is searched for a base object,
   * which is then again stored.
   *
   * @param gridConfig The desired grid configuration.
   * @param geometryConfig The desired geometry configuration.
   * @param adaptConfig The desired adaptivity configuration.
   * @param regularizationConfig The desired regularization configuration.
   * @param densityEstimationConfig The desired desnity estimation configuration.
   * @return A offline object that matches the configuration.
   */
  DBMatOfflinePermutable* getPermutedObject(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

 protected:
  std::shared_ptr<DBMatObjectStore> store;
  bool hasDataBase;
  std::string dbFilePath;
};
}  // namespace datadriven
}  // namespace sgpp
