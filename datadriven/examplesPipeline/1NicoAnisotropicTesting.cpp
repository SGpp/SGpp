#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

int main(int argc, char** argv) {
  const std::string path = [argc, &argv]() {
    if (argc != 2) {
      std::cout << "No or bad path given, aborting\n";
      exit(1);
      return std::string{};
    } else {
      return std::string{argv[1]};
    }
  }();

  const sgpp::datadriven::DataMiningConfigParser parser(path);
  sgpp::datadriven::FitterConfigurationDensityEstimation config{};
  config.readParams(parser);

  size_t dim = config.getGridConfig().dim_;
  size_t level = config.getGridConfig().level_;
  std::vector<size_t> levelvec = config.getGridConfig().levelVector_;

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  grid->getGenerator().anisotropicFull(levelvec);

  std::cout << "all done\n";
  return 1;
}
