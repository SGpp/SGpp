#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>

namespace sgpp {
namespace datadriven {

class DBMatBaseObjectStore {
 public:
  DBMatBaseObjectStore();

  explicit DBMatBaseObjectStore(const std::string& fileName);

  /**
   * Returns a offline object for the given configuration.
   * If no base object is currently available, the method first tries to obtain
   * a suitable decomposition from the database, if a database file is specified.
   * If no database is specified or the database contains no suitable base object,
   * a new base object is built and decomposed.
   * Note that base object allways don't have level vector elements equal to 1.
   * Hence, the dimension blow up method may be apllied to a base object
   * constructed from the database.
   * The base object then gets transformed into the desired offline object,
   * by applying permutation and dimension blow-up.
   * @param gridConfig the desired component grid configuration
   * @param adaptivityConfig the desired adaptivity configuration
   * @param regularizationConfig the desired regularization configuration
   * @param densityEstimationConfig the desired desityestimation config
   */
  DBMatOfflinePermutable* getOfflineObject(sgpp::base::GeneralGridConfiguration& gridConfig,
                                           sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                                           RegularizationConfiguration& regularizationConfig,
                                           DensityEstimationConfiguration& densityEstimationConfig);

 protected:
  class ObjectContainer {
   public:
    explicit ObjectContainer(const sgpp::base::GeneralGridConfiguration& gridConfig,
                             const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                             const RegularizationConfiguration& regularizationConfig,
                             const DensityEstimationConfiguration& densityEstimationConfig,
                             std::unique_ptr<const DBMatOfflinePermutable> offlineObject);

    const DBMatOfflinePermutable& getOfflineObject() const;

    const sgpp::base::GeneralGridConfiguration& getGridConfig() const;

    bool configMatches(const sgpp::base::GeneralGridConfiguration& gridConfig,
                       const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                       const RegularizationConfiguration& regularizationConfig,
                       const DensityEstimationConfiguration& densityEstimationConfig);

   private:
    sgpp::base::GeneralGridConfiguration gridConfig;
    sgpp::base::AdaptivityConfiguration adaptivityConfig;
    RegularizationConfiguration regularizationConfig;
    DensityEstimationConfiguration densityEstimationConfig;
    std::unique_ptr<const DBMatOfflinePermutable> offlineObject;
  };

  std::vector<ObjectContainer> objects;
  std::string dbFilePath;
  bool hasDatabase;
  int getBaseObjectContainerIndex(sgpp::base::GeneralGridConfiguration& gridConfig,
                                  sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                                  RegularizationConfiguration& regularizationConfig,
                                  DensityEstimationConfiguration& densityEstimationConfig);
  const ObjectContainer& getBaseObjectContainer(int index) const;
  void putBaseObject(sgpp::base::GeneralGridConfiguration& gridConfig,
                     sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                     RegularizationConfiguration& regularizationConfig,
                     DensityEstimationConfiguration& densityEstimationConfig,
                     std::unique_ptr<DBMatOfflinePermutable> object);
};

}  // namespace datadriven
}  // namespace sgpp