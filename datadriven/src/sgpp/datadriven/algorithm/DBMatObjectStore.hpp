#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>

namespace sgpp {
namespace datadriven {

class DBMatObjectStore {
 public:
  DBMatObjectStore();

  explicit DBMatObjectStore(const std::string& fileName);

  void putObject(const sgpp::base::GeneralGridConfiguration& gridConfig,
                 const sgpp::datadriven::GeometryConfiguration& geometryConfig,
                 const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                 const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                 const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
                 const DBMatOffline* object);

  const DBMatOfflinePermutable* getBaseObject(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration& geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      sgpp::base::GeneralGridConfiguration& baseGridConfig);

  const DBMatOffline* getObject(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration& geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

 protected:
  class ObjectContainer {
   public:
    explicit ObjectContainer(
        const sgpp::base::GeneralGridConfiguration& gridConfig,
        const sgpp::datadriven::GeometryConfiguration& geometryConfig,
        const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
        const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
        const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
        std::unique_ptr<const DBMatOffline> offlineObject);

    const DBMatOffline& getOfflineObject() const;

    const sgpp::base::GeneralGridConfiguration& getGridConfig() const;

    bool configMatches(const sgpp::base::GeneralGridConfiguration& gridConfig,
                       const sgpp::datadriven::GeometryConfiguration& geometryConfig,
                       const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                       const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                       const sgpp::datadriven::DensityEstimationConfiguration&
                           densityEstimationConfig, bool searchBase);

   private:
    sgpp::base::GeneralGridConfiguration gridConfig;
    sgpp::datadriven::GeometryConfiguration geometryConfig;
    sgpp::base::AdaptivityConfiguration adaptivityConfig;
    RegularizationConfiguration regularizationConfig;
    DensityEstimationConfiguration densityEstimationConfig;
    std::unique_ptr<const DBMatOffline> offlineObject;
  };

  std::vector<ObjectContainer> objects;
  std::string dbFilePath;
  bool hasDatabase;
  int getObjectContainerIndex(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration& geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
      bool searchBase = false);
  const ObjectContainer& getObjectContainer(int index) const;
};

}  // namespace datadriven
}  // namespace sgpp