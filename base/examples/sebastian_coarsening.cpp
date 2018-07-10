#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>

#include <iostream>
#include <algorithm>

using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::HashCoarsening;
using sgpp::base::SurplusCoarseningFunctor;


void printLevel(GridPoint& gp) {
  for(size_t i = 0; i < 3; i++) {
    std::cout << gp.getLevel(i) << " ";
  }
  std::cout << "| i=" << gp.isInnerPoint() << " l=" << gp.isLeaf();
  if(gp.isInnerPoint() && gp.isLeaf()) {
    std::cout << " < coarsable";
  }
  if(gp.isLeaf() && !gp.isInnerPoint()) {
    std::cout << " < what?";
  }
  std::cout << std::endl;
}

int main() {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(Grid::createLinearBoundaryGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  std::cout << "Number of GP: " << gridStorage.getSize() << std::endl;

  /*
  for(size_t i = 0; i < gridStorage.getSize(); i++) {
    std::cout << i << ": ";
    printLevel(gridStorage.getPoint(i));
  }
  */

  std::vector<size_t> toBeRemoved;
  DataVector alpha(gridStorage.getSize());
  alpha.setAll(1.0);
  alpha[161] = 0.1;  // should not be coarsed i=0 l=0
  alpha[162] = 0.1;  // should not be coarsed i=1 l=0
  alpha[203] = 0.1;  toBeRemoved.push_back(gridStorage.getPoint(203).getHash()); // should be coarsed i=1 l=1
  alpha[204] = 0.1;  toBeRemoved.push_back(gridStorage.getPoint(204).getHash()); // should be coarsed i=1 l=1
  alpha[208] = 0.1;  toBeRemoved.push_back(gridStorage.getPoint(208).getHash()); // should be coarsed i=1 l=1

  for(size_t i = 0; i < toBeRemoved.size(); i++) {
    std::cout << toBeRemoved[i] << std::endl;
  }

  std::vector<size_t> before;
  for(auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    before.push_back(it->first->getHash());
  }

  std::cout << "Size before: " << gridStorage.getSize() << std::endl;
  HashCoarsening coarsen;
  SurplusCoarseningFunctor functor(alpha, 3, 0.5);
  coarsen.free_coarsen(gridStorage, functor, alpha);
  std::cout << "Size after: " << gridStorage.getSize() << std::endl;

  std::vector<size_t> after;
  for(auto it = gridStorage.begin(); it != gridStorage.end(); it++) {
    after.push_back(it->first->getHash());
  }

  for(size_t i = 0; i < before.size(); i++) {
    if(std::find(after.begin(), after.end(), before.at(i)) == after.end()) {
      if(std::find(toBeRemoved.begin(), toBeRemoved.end(), before.at(i)) != toBeRemoved.end()) {
        std::cout << "Point was removed" << std::endl;
      } else {
        std::cout << "Whaaaa" << std::endl;
      }
    }
  }

}
