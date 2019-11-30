// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/tools/vpTree/VpTree.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>
#include <iostream>

using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::datadriven::VpTree;

BOOST_AUTO_TEST_SUITE(test_vpTree)

BOOST_AUTO_TEST_CASE(nearestNeighbors) {
  // Test 1. Deliver the correct number of nearest neighbors
  double points[] = {0.5, 0, -0.6, 0, 0, 0.3, 0, -0.2, 0.1, 0.2, 0, 0};
  double origin[] = {0, 0};

  DataVector target(origin, 2);

  DataMatrix inputMatrix(points, 6, 2);

  std::unique_ptr<VpTree> tree = std::make_unique<VpTree>(inputMatrix);

  auto nearestNeighborsHeap = tree->getNearestNeighbors(target, 3);

  std::vector<size_t> nearestNeighbors;

  nearestNeighbors.clear();
  while (!nearestNeighborsHeap.empty()) {
        nearestNeighbors.push_back(nearestNeighborsHeap.top().index);
        nearestNeighborsHeap.pop();
  }

  std::reverse(nearestNeighbors.begin(), nearestNeighbors.end());

  std::cout << "Number of nearest Neighbors "
  << std::to_string(nearestNeighbors.size()) << std::endl;

  BOOST_CHECK_EQUAL(nearestNeighbors.size(), 3);

  // Test 2. Deliver the correct nearest neighbors of the following points
  DataVector neighbor1(2);
  double point1[] = {0, -0.2};
  DataVector trueNeighbor1(point1, 2);
  tree->getStoredItems().getRow(nearestNeighbors.at(0), neighbor1);
  std::cout << "Neighbors 1 " << neighbor1.toString() << std::endl;
  std::cout << "True Neighbor 1 " << trueNeighbor1.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor1.get(0), trueNeighbor1.get(0));
  BOOST_CHECK_EQUAL(neighbor1.get(1), trueNeighbor1.get(1));

  DataVector neighbor2(2);
  double point2[] = {0.1, 0.2};
  DataVector trueNeighbor2(point2, 2);
  tree->getStoredItems().getRow(nearestNeighbors.at(1), neighbor2);
  std::cout << "Neighbors 2 " << neighbor2.toString() << std::endl;
  std::cout << "True Neighbor 2 " << trueNeighbor2.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor2.get(0), trueNeighbor2.get(0));
  BOOST_CHECK_EQUAL(neighbor2.get(1), trueNeighbor2.get(1));

  DataVector neighbor3(2);
  double point3[] = {0, 0.3};
  DataVector trueNeighbor3(point3, 2);
  tree->getStoredItems().getRow(nearestNeighbors.at(2), neighbor3);
  std::cout << "Neighbors 3 " << neighbor3.toString() << std::endl;
  std::cout << "True Neighbor 3 " << trueNeighbor3.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor3.get(0), trueNeighbor3.get(0));
  BOOST_CHECK_EQUAL(neighbor3.get(1), trueNeighbor3.get(1));

  // Test 3. Deliver the same nearest neighbors of the same points after and online update if
  // the given points are no new nearest neighbors
  double updatePoint[] = {1.0, 0};
  DataMatrix updateMatrix(updatePoint, 1, 2);
  tree->update(updateMatrix);

  nearestNeighborsHeap = tree->getNearestNeighbors(target, 3);

  nearestNeighbors.clear();
  while (!nearestNeighborsHeap.empty()) {
        nearestNeighbors.push_back(nearestNeighborsHeap.top().index);
        nearestNeighborsHeap.pop();
  }

  std::reverse(nearestNeighbors.begin(), nearestNeighbors.end());

  tree->getStoredItems().getRow(nearestNeighbors.at(0), neighbor1);
  std::cout << "Neighbors 1 " << neighbor1.toString() << std::endl;
  std::cout << "True Neighbor 1 " << trueNeighbor1.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor1.get(0), trueNeighbor1.get(0));
  BOOST_CHECK_EQUAL(neighbor1.get(1), trueNeighbor1.get(1));

  tree->getStoredItems().getRow(nearestNeighbors.at(1), neighbor2);
  std::cout << "Neighbors 2 " << neighbor2.toString() << std::endl;
  std::cout << "True Neighbor 2 " << trueNeighbor2.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor2.get(0), trueNeighbor2.get(0));
  BOOST_CHECK_EQUAL(neighbor2.get(1), trueNeighbor2.get(1));

  tree->getStoredItems().getRow(nearestNeighbors.at(2), neighbor3);
  std::cout << "Neighbors 3 " << neighbor3.toString() << std::endl;
  std::cout << "True Neighbor 3 " << trueNeighbor3.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor3.get(0), trueNeighbor3.get(0));
  BOOST_CHECK_EQUAL(neighbor3.get(1), trueNeighbor3.get(1));


  // Test 4. Deliver the new nearest neighbors of the same points after and online update if
  // the given update points are indeed new nearest neighbors
  updateMatrix.set(0, 0, 0.1);
  updateMatrix.set(0, 1, 0.001);
  tree->update(updateMatrix);
  std::cout << tree->getStoredItems().toString() << std::endl;
  nearestNeighborsHeap = tree->getNearestNeighbors(target, 3);

  nearestNeighbors.clear();
  while (!nearestNeighborsHeap.empty()) {
        nearestNeighbors.push_back(nearestNeighborsHeap.top().index);
        nearestNeighborsHeap.pop();
  }

  std::reverse(nearestNeighbors.begin(), nearestNeighbors.end());

  trueNeighbor3.copyFrom(trueNeighbor2);
  trueNeighbor2.copyFrom(trueNeighbor1);
  trueNeighbor1.set(0, 0.1);
  trueNeighbor1.set(1, 0.001);
  tree->getStoredItems().getRow(nearestNeighbors.at(0), neighbor1);
  std::cout << "Neighbors 1 " << neighbor1.toString() << std::endl;
  std::cout << "True Neighbor 1 " << trueNeighbor1.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor1.get(0), trueNeighbor1.get(0));
  BOOST_CHECK_EQUAL(neighbor1.get(1), trueNeighbor1.get(1));

  tree->getStoredItems().getRow(nearestNeighbors.at(1), neighbor2);
  std::cout << "Neighbors 2 " << neighbor2.toString() << std::endl;
  std::cout << "True Neighbor 2 " << trueNeighbor2.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor2.get(0), trueNeighbor2.get(0));
  BOOST_CHECK_EQUAL(neighbor2.get(1), trueNeighbor2.get(1));

  tree->getStoredItems().getRow(nearestNeighbors.at(2), neighbor3);
  std::cout << "Neighbors 3 " << neighbor3.toString() << std::endl;
  std::cout << "True Neighbor 3 " << trueNeighbor3.toString() << std::endl;
  BOOST_CHECK_EQUAL(neighbor3.get(0), trueNeighbor3.get(0));
  BOOST_CHECK_EQUAL(neighbor3.get(1), trueNeighbor3.get(1));
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
