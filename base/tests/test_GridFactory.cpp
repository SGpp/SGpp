// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <algorithm>
#include <vector>
#include <string>

using sgpp::base::BoundingBox;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::base::BoundingBox1D;
using sgpp::base::generation_exception;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::OperationEval;
using sgpp::base::OperationMultipleEval;
using sgpp::base::Stretching;
using sgpp::base::Stretching1D;
using sgpp::base::SurplusRefinementFunctor;

BOOST_AUTO_TEST_SUITE(TestGridFactory)

BOOST_AUTO_TEST_CASE(testCreation) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createLinearGrid(2));
  BOOST_CHECK(factory != nullptr);
}

BOOST_AUTO_TEST_CASE(testSerializationLinear) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createLinearGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());
}

BOOST_AUTO_TEST_CASE(testSerializationModLinear) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createModLinearGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());
}

BOOST_AUTO_TEST_CASE(testSerializationLinearTruncatedBoundary) {
  // Uses Linear grid for tests
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());
}

BOOST_AUTO_TEST_CASE(testSerializationLinearBoundary) {
  // Uses Linear grid for tests
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());
}

BOOST_AUTO_TEST_CASE(testSerializationPrewavelet) {
  // Uses Linear grid for tests
  std::unique_ptr<Grid> factory(Grid::createPrewaveletGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());
}

BOOST_AUTO_TEST_CASE(testSerializationLinearBoundingBox) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createLinearGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  {
    BoundingBox& boundingBox = factory->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox.setBoundary(0, tempBound);
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  {
    BoundingBox& boundingBox = newfac->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    BOOST_CHECK(0.0 == tempBound.leftBoundary);
    BOOST_CHECK(100.0 == tempBound.rightBoundary);
    BOOST_CHECK(false == tempBound.bDirichletLeft);
    BOOST_CHECK(false == tempBound.bDirichletRight);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationModLinearBoundingBox) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createModLinearGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  {
    BoundingBox& boundingBox = factory->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox.setBoundary(0, tempBound);
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  {
    BoundingBox& boundingBox = newfac->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    BOOST_CHECK(0.0 == tempBound.leftBoundary);
    BOOST_CHECK(100.0 == tempBound.rightBoundary);
    BOOST_CHECK(false == tempBound.bDirichletLeft);
    BOOST_CHECK(false == tempBound.bDirichletRight);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationLinearTruncatedBoundaryBoundingBox) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  {
    BoundingBox& boundingBox = factory->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox.setBoundary(0, tempBound);
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  {
    BoundingBox& boundingBox = newfac->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    BOOST_CHECK(0.0 == tempBound.leftBoundary);
    BOOST_CHECK(100.0 == tempBound.rightBoundary);
    BOOST_CHECK(false == tempBound.bDirichletLeft);
    BOOST_CHECK(false == tempBound.bDirichletRight);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationLinearBoundaryBoundingBox) {
  // Uses Linear grid for tests

  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  {
    BoundingBox& boundingBox = factory->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox.setBoundary(0, tempBound);
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  {
    BoundingBox& boundingBox = newfac->getBoundingBox();
    BoundingBox1D tempBound = boundingBox.getBoundary(0);
    BOOST_CHECK(0.0 == tempBound.leftBoundary);
    BOOST_CHECK(100.0 == tempBound.rightBoundary);
    BOOST_CHECK(false == tempBound.bDirichletLeft);
    BOOST_CHECK(false == tempBound.bDirichletRight);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationLinearWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  std::unique_ptr<Grid> factory(Grid::createLinearGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  for (size_t i = 0; i < factory->getSize(); ++i) {
    srcLeaf.push_back(factory->getStorage().getPoint(i).isLeaf());
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  for (size_t i = 0; i < factory->getSize(); ++i) {
    BOOST_CHECK(newfac->getStorage().getPoint(i).isLeaf() == srcLeaf[i]);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationModLinearWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  std::unique_ptr<Grid> factory(Grid::createModLinearGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  for (size_t i = 0; i < factory->getSize(); ++i) {
    srcLeaf.push_back(factory->getStorage().getPoint(i).isLeaf());
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  for (size_t i = 0; i < factory->getSize(); ++i) {
    BOOST_CHECK(newfac->getStorage().getPoint(i).isLeaf() == srcLeaf[i]);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationLinearTruncatedBoundaryWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  for (size_t i = 0; i < factory->getSize(); ++i) {
    srcLeaf.push_back(factory->getStorage().getPoint(i).isLeaf());
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  for (size_t i = 0; i < factory->getSize(); ++i) {
    BOOST_CHECK(newfac->getStorage().getPoint(i).isLeaf() == srcLeaf[i]);
  }
}

BOOST_AUTO_TEST_CASE(testSerializationLinearBoundaryWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  BOOST_CHECK(factory != nullptr);

  factory->getGenerator().regular(3);

  for (size_t i = 0; i < factory->getSize(); ++i) {
    srcLeaf.push_back(factory->getStorage().getPoint(i).isLeaf());
  }

  std::string str = factory->serialize();
  BOOST_CHECK(str.size() > 0);

  std::unique_ptr<Grid> newfac(Grid::unserialize(str));
  BOOST_CHECK(newfac != nullptr);

  BOOST_CHECK(factory->getSize() == newfac->getSize());

  for (size_t i = 0; i < factory->getSize(); ++i) {
    BOOST_CHECK(newfac->getStorage().getPoint(i).isLeaf() == srcLeaf[i]);
  }
}

// end test suite TestGridFactory
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  std::unique_ptr<Grid> factory(Grid::createLinearGrid(2));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();

  BOOST_CHECK(storage.getSize() == 0);
  gen.regular(3);
  BOOST_CHECK(storage.getSize() == 17);

  // This should fail
  BOOST_CHECK_THROW(gen.regular(3), generation_exception);
}

BOOST_AUTO_TEST_CASE(testRefinement) {
  std::unique_ptr<Grid> factory(Grid::createLinearGrid(2));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  BOOST_CHECK(storage.getSize() == 1);
  DataVector alpha(1);
  alpha[0] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);
  BOOST_CHECK(storage.getSize() == 5);
}

BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  std::unique_ptr<Grid> factory(Grid::createLinearGrid(1));
  GridGenerator& gen = factory->getGenerator();
  gen.regular(2);

  DataVector alpha(factory->getSize(), 0.0);
  DataMatrix p(1, 1);
  DataVector beta(1);

  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  std::unique_ptr<OperationMultipleEval> opb(
      sgpp::op_factory::createOperationMultipleEval(*factory, p));
  opb->multTranspose(beta, alpha);

  BOOST_CHECK_CLOSE(alpha[0], 0.5, 0.0);
  BOOST_CHECK_CLOSE(alpha[1], 1.0, 0.0);
  BOOST_CHECK_CLOSE(alpha[2], 0.0, 0.0);

  alpha.setAll(0.0);
  alpha[0] = 1.0;

  p.set(0, 0, 0.25);

  beta[0] = 0.0;

  opb->mult(alpha, beta);
  BOOST_CHECK_CLOSE(beta[0], 0.5, 0.0);
}

BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  std::unique_ptr<Grid> factory(Grid::createLinearGrid(1));
  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  DataVector alpha(factory->getSize(), 1.0);

  DataVector p(1);
  p.setAll(0.25);

  std::unique_ptr<OperationEval> eval(
      sgpp::op_factory::createOperationEval(*factory));
  BOOST_CHECK_CLOSE(eval->eval(alpha, p), 0.5, 0.0);
}

// end test suite TestLinearGrid
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearBoundaryGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();

  BOOST_CHECK(storage.getSize() == 0);
  gen.regular(3);
  BOOST_CHECK(storage.getSize() == 49);

  // This should fail
  BOOST_CHECK_THROW(gen.regular(3), generation_exception);
}

BOOST_AUTO_TEST_CASE(testRefinement2d) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  BOOST_CHECK(storage.getSize() == 9);
  DataVector alpha(9);
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 0.0;
  alpha[3] = 0.0;
  alpha[4] = 0.0;
  alpha[5] = 0.0;
  alpha[6] = 0.0;
  alpha[7] = 0.0;
  alpha[8] = 1.0;

  SurplusRefinementFunctor func(alpha);

  gen.refine(func);
  BOOST_CHECK(storage.getSize() == 21);
}

BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(3));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  BOOST_CHECK(storage.getSize() == 27);
  DataVector alpha(27);
  alpha.setAll(0.0);

  alpha[26] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);
  BOOST_CHECK(storage.getSize() == 81);
}

BOOST_AUTO_TEST_CASE(testRefinement3d) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(1));
  GridGenerator& gen = factory->getGenerator();
  gen.regular(2);

  DataVector alpha(factory->getSize(), 0.0);
  DataMatrix p(1, 1);
  DataVector beta(1);

  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  std::unique_ptr<OperationMultipleEval> opb(
      sgpp::op_factory::createOperationMultipleEval(*factory, p));
  opb->multTranspose(beta, alpha);

  BOOST_CHECK_CLOSE(alpha[0], 0.75, 0.0);
  BOOST_CHECK_CLOSE(alpha[1], 0.25, 0.0);
  BOOST_CHECK_CLOSE(alpha[2], 0.5, 0.0);
  BOOST_CHECK_CLOSE(alpha[3], 1.0, 0.0);
  BOOST_CHECK_CLOSE(alpha[4], 0.0, 0.0);

  alpha.setAll(0.0);
  alpha[2] = 1.0;

  p.set(0, 0, 0.25);

  beta[0] = 0.0;

  opb->mult(alpha, beta);
  BOOST_CHECK_CLOSE(beta[0], 0.5, 0.0);
}

BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(1));
  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  DataVector alpha(factory->getSize(), 1.0);

  DataVector p(1, 0.25);

  std::unique_ptr<OperationEval> eval(
      sgpp::op_factory::createOperationEval(*factory));

  BOOST_CHECK_CLOSE(eval->eval(alpha, p), 1.5, 0.0);
}

// end test suite TestLinearBoundaryGrid
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearL0BoundaryGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();

  BOOST_CHECK(storage.getSize() == 0);
  gen.regular(3);
  BOOST_CHECK(storage.getSize() == 37);

  // This should fail
  BOOST_CHECK_THROW(gen.regular(3), generation_exception);
}

BOOST_AUTO_TEST_CASE(testRefinement2d_one) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(0);

  BOOST_CHECK(storage.getSize() == 4);

  DataVector alpha(4);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);
  BOOST_CHECK(storage.getSize() == 8);
}

BOOST_AUTO_TEST_CASE(testRefinement2d_two) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(0);

  DataVector alpha(4);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);

  DataVector alpha2(8);

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2(alpha2);

  gen.refine(func2);
  BOOST_CHECK(storage.getSize() == 13);
}

BOOST_AUTO_TEST_CASE(testRefinement2d_three) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(0);

  DataVector alpha(4);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);

  DataVector alpha2(8);

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2(alpha2);

  gen.refine(func2);

  DataVector alpha3(13);

  alpha3[11] = 1.0;
  SurplusRefinementFunctor func3(alpha3);

  gen.refine(func3);
  BOOST_CHECK(storage.getSize() == 18);
}

BOOST_AUTO_TEST_CASE(testRefinement2d_four) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(0);

  DataVector alpha(4);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);

  DataVector alpha2(8);

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2(alpha2);

  gen.refine(func2);

  DataVector alpha3(13);

  alpha3[11] = 1.0;
  SurplusRefinementFunctor func3(alpha3);

  gen.refine(func3);

  DataVector alpha4(18);

  alpha4[12] = 1.0;
  SurplusRefinementFunctor func4(alpha4);

  gen.refine(func4);
  BOOST_CHECK(storage.getSize() == 25);
}

BOOST_AUTO_TEST_CASE(testRefinement2d_five) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(0);

  DataVector alpha(4);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);

  DataVector alpha2(8);

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2(alpha2);

  gen.refine(func2);

  DataVector alpha3(13);

  alpha3[11] = 1.0;
  SurplusRefinementFunctor func3(alpha3);

  gen.refine(func3);

  DataVector alpha4(18);

  alpha4[12] = 1.0;
  SurplusRefinementFunctor func4(alpha4);

  gen.refine(func4);

  DataVector alpha5(25);

  alpha5[23] = 1.0;
  SurplusRefinementFunctor func5(alpha5);

  gen.refine(func5);

  BOOST_CHECK(storage.getSize() == 29);
}

BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridGenerator& gen = factory->getGenerator();
  gen.regular(2);
  // point 0: l = (0,0), i = (0,0)
  // point 1: l = (0,0), i = (1,0)
  // point 2: l = (1,0), i = (1,0)
  // point 3: l = (2,0), i = (1,0)
  // point 4: l = (2,0), i = (3,0)

  DataVector alpha(factory->getSize(), 0.0);
  DataMatrix p(1, 2);
  DataVector beta(1);

  p.set(0, 0, 0.25);
  p.set(0, 1, 0.25);
  beta[0] = 1.0;

  std::unique_ptr<OperationMultipleEval> opb(
      sgpp::op_factory::createOperationMultipleEval(*factory, p));
  opb->multTranspose(beta, alpha);

  // should be phi_{(0,0),(0,0)}(p) = (3/4)^2 = 9/16
  BOOST_CHECK_CLOSE(alpha[0], 0.5625, 0.0);
  // should be phi_{(0,0),(1,0)}(p) = (1/4)*(3/4) = 3/16
  BOOST_CHECK_CLOSE(alpha[1], 0.1875, 0.0);
  // should be phi_{(1,0),(1,0)}(p) = (1/2)*(3/4) = 3/8
  BOOST_CHECK_CLOSE(alpha[2], 0.375, 0.0);
  // should be phi_{(2,0),(1,0)}(p) = 1*3/4 = 3/4
  BOOST_CHECK_CLOSE(alpha[3], 0.75, 0.0);
  // should be phi_{(2,0),(3,0)}(p) = 0*3/4 = 0
  BOOST_CHECK_CLOSE(alpha[4], 0.0, 0.0);

  alpha.setAll(0.0);
  alpha[2] = 1.0;

  beta[0] = 0.0;

  opb->mult(alpha, beta);

  // rationale behind 0.375: function no. 2 should have level
  // l = (1,0) and index i = (1,0)
  // ==>  function evaluated at (0.25, 0.25) should be 1/2 * 3/4 = 3/8
  BOOST_CHECK_CLOSE(beta[0], 0.375, 0.0);
}

BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  std::unique_ptr<Grid> factory(Grid::createLinearBoundaryGrid(2, 0));
  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  DataVector alpha(factory->getSize(), 1.0);

  DataVector p(2, 0.25);

  std::unique_ptr<OperationEval> eval(
      sgpp::op_factory::createOperationEval(*factory));

  // rationale behind 2.0: four corner functions sum up to 1,
  // two of the four edge functions are 1/8, the other two are 3/8
  // ==>  1 + 2 * 1/8 + 2 * 3/8 = 2
  BOOST_CHECK_CLOSE(eval->eval(alpha, p), 2.0, 0.0);
}

// end test suite TestLinearL0BoundaryGrid
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearStretchedBoundaryGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;

  BoundingBox1D dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;

  std::vector<BoundingBox1D> dimBoundaryVector(2);
  dimBoundaryVector[0] = dimBound;
  dimBoundaryVector[1] = dimBound;

  std::vector<Stretching1D> str1dvector(2);
  str1dvector[0] = str1d;
  str1dvector[1] = str1d;
  Stretching stretch(dimBoundaryVector, str1dvector);

  std::unique_ptr<Grid> factory(Grid::createLinearStretchedBoundaryGrid(2));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();

  BOOST_CHECK(storage.getSize() == 0);
  gen.regular(3);
  BOOST_CHECK(storage.getSize() == 49);

  // This should fail
  BOOST_CHECK_THROW(gen.regular(3), generation_exception);
}

BOOST_AUTO_TEST_CASE(testRefinement2d) {
  std::unique_ptr<Grid> factory(Grid::createLinearStretchedBoundaryGrid(2));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  BOOST_CHECK(storage.getSize() == 9);
  DataVector alpha(9);
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 0.0;
  alpha[3] = 0.0;
  alpha[4] = 0.0;
  alpha[5] = 0.0;
  alpha[6] = 0.0;
  alpha[7] = 0.0;
  alpha[8] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);
  BOOST_CHECK(storage.getSize() == 21);
}

BOOST_AUTO_TEST_CASE(testRefinement3d) {
  std::unique_ptr<Grid> factory(Grid::createLinearStretchedBoundaryGrid(3));
  GridStorage& storage = factory->getStorage();

  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  BOOST_CHECK(storage.getSize() == 27);
  DataVector alpha(27);

  alpha[26] = 1.0;
  SurplusRefinementFunctor func(alpha);

  gen.refine(func);
  BOOST_CHECK(storage.getSize() == 81);
}

BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;
  BoundingBox1D dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;
  Stretching stretch({dimBound}, {str1d});

  std::unique_ptr<Grid> factory(Grid::createLinearStretchedBoundaryGrid(1));
  factory->getStorage().setStretching(stretch);
  GridGenerator& gen = factory->getGenerator();
  gen.regular(2);

  DataVector alpha(factory->getSize(), 0.0);

  DataMatrix p(1, 1);
  DataVector beta(1);

  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  std::unique_ptr<OperationMultipleEval> opb(
      sgpp::op_factory::createOperationMultipleEval(*factory, p));
  opb->multTranspose(beta, alpha);

  BOOST_CHECK_CLOSE(alpha[0], 1.038461538, 1e-6);
  BOOST_CHECK_CLOSE(alpha[1], -0.038461538461538464, 1e-6);
  BOOST_CHECK_CLOSE(alpha[2], -0.18237143795284394, 1e-6);
  BOOST_CHECK_CLOSE(alpha[3], -0.53513915, 1e-6);
  BOOST_CHECK_CLOSE(alpha[4], 0.0, 1e-6);

  alpha.setAll(0.0);
  alpha[2] = 1.0;

  p.set(0, 0, 0.25);

  beta[0] = 0.0;

  opb->mult(alpha, beta);
  BOOST_CHECK_CLOSE(beta[0], -0.182371437, 1e-6);
}

BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  /*
  from pysgpp import Grid, DataVector, Stretching, Stretching1D, BoundingBox1D

  str1d = Stretching1D()
  str1d.type='log'
  str1d.x_0=1
  str1d.xsi=10
  dimBound = BoundingBox1D()
  dimBound.leftBoundary=0.5
  dimBound.rightBoundary=7
  stretch=Stretching(1,dimBound,str1d)

  factory = Grid.createLinearStretchedBoundaryGrid(1)
  factory.getStorage().setStretching(stretch)
  gen = factory.getGenerator()
  gen.regular(1)

  alpha = DataVector(factory.getSize(), 1.0)

  p = DataVector(1, 0.25)

  eval = createOperationEval(factory)

  self.failUnlessAlmostEqual(eval.eval(alpha, p), 0.8176285620)
  */
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;
  BoundingBox1D dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;
  Stretching stretch({dimBound}, {str1d});

  std::unique_ptr<Grid> factory(Grid::createLinearStretchedBoundaryGrid(1));
  factory->getStorage().setStretching(stretch);
  GridGenerator& gen = factory->getGenerator();
  gen.regular(1);

  DataVector alpha(factory->getSize(), 1.0);

  DataVector p(1, 0.25);

  std::unique_ptr<OperationEval> eval(
      sgpp::op_factory::createOperationEval(*factory));

  BOOST_CHECK_CLOSE(eval->eval(alpha, p), 0.8176285620, 1e-8);
}
// end test suite TestLinearStretchedBoundaryGrid
BOOST_AUTO_TEST_SUITE_END()
