#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflinePermutable.hpp>

namespace sgpp {
namespace datadriven {
class DBMatPermutationFactory {
 public:
  DBMatPermutationFactory();
  explicit DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore> store);
  explicit DBMatPermutationFactory(std::shared_ptr<DBMatObjectStore>  store, const std::string& dbFilePath);

  DBMatOfflinePermutable* getPermutedObject(
      const sgpp::base::GeneralGridConfiguration& gridConfig,
      const sgpp::datadriven::GeometryConfiguration geometryConfig,
      const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
      const sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      const sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig);

 protected:
  std::string dbFilePath;
  bool hasDataBase;
  std::shared_ptr<DBMatObjectStore> store;
};
}  // namespace datadriven
}  // namespace sgpp