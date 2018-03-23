//
// Created by Eric on 20.03.2018.
//

#ifndef CLION_BOCONFIG_HPP
#define CLION_BOCONFIG_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <random>

namespace sgpp {
namespace datadriven {

class BOConfig {
public:
  BOConfig() = default;
  BOConfig(std::vector<int>* discOptions, std::vector<int>* catOptions, size_t nCont);

  bool nextDisc();

  void calcDiscDistance(BOConfig& other);

  double getTotalDistance(const base::DataVector& input);

  size_t getContSize();

  void setCont(const base::DataVector& input);

  double getCont(int idx);

  int getDisc(int idx);

  int getCat(int idx);

  void setScore(double input);

  double getScore();

  double getDistance(BOConfig& other);

  void randomize(std::mt19937& generator);


private:
  base::DataVector cont;
  std::vector<int> disc;
  std::vector<int> cat;
  std::vector<int>* discOptions = nullptr;
  std::vector<int>* catOptions = nullptr;
  double score = 0;
  double discDistance = 0;
};

}
}

#endif //CLION_BOCONFIG_HPP
