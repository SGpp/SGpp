// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>

#include <string>

using sgpp::base::DataVector;
using sgpp::base::HashGenerator;
using sgpp::base::HashGridIndex;
using sgpp::base::HashGridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::HashRefinementBoundaries;
using sgpp::base::SurplusRefinementFunctor;

BOOST_AUTO_TEST_SUITE(TestHashGridStorage)

BOOST_AUTO_TEST_CASE(testCreateDestroy) {
  HashGridIndex i(1);
  HashGridStorage s(1);

  i.set(0, 1, 1);
  HashGridIndex* i2 = s.create(i);

  HashGridIndex::level_type l, l2;
  HashGridIndex::index_type ind, ind2;

  i.get(0, l, ind);
  i2->get(0, l2, ind2);

  BOOST_CHECK_EQUAL(ind, ind2);
  BOOST_CHECK_EQUAL(l, l2);

  s.destroy(i2);
}

BOOST_AUTO_TEST_CASE(testSerialize) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regular(s, 2);

  std::string str = s.serialize();

  BOOST_TEST_MESSAGE(str);
  BOOST_CHECK_NE(str.length(), 0U);

  HashGridStorage s2(str);

  BOOST_CHECK_EQUAL(s.getSize(), s2.getSize());
}

BOOST_AUTO_TEST_CASE(testSerializeWithLeaf) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regular(s, 2);

  bool* srcLeaf = new bool[s.getSize()];

  for (unsigned int i = 0; i < s.getSize(); ++i) {
    srcLeaf[i] = s.get(i)->isLeaf();
  }

  std::string str = s.serialize();

  BOOST_CHECK_NE(str.length(), 0U);

  HashGridStorage s2(str);

  BOOST_CHECK_EQUAL(s.getSize(), s2.getSize());

  for (unsigned int i = 0; i < s.getSize(); ++i) {
    BOOST_CHECK_EQUAL(s2.get(i)->isLeaf(), srcLeaf[i]);
  }

  delete[] srcLeaf;
}

BOOST_AUTO_TEST_CASE(testInsert) {
  HashGridIndex i(1);
  HashGridStorage s(1);

  i.set(0, 1, 1);
  size_t i2 = s.insert(i);

  BOOST_CHECK_EQUAL(i2, 0U);
  BOOST_CHECK_EQUAL(s.getSize(), 1U);
}

BOOST_AUTO_TEST_CASE(testChilds) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regular(s, 2);

  HashGridIndex i(1);

  i.set(0, 1, 1);
  s.left_child(&i, 0);

  HashGridIndex::level_type l, l2 = 2;
  HashGridIndex::index_type ind, ind2 = 1;

  i.get(0, l, ind);
  BOOST_CHECK_EQUAL(l, l2);
  BOOST_CHECK_EQUAL(ind, ind2);

  i.set(0, 1, 1);
  s.right_child(&i, 0);
  l2 = 2;
  ind2 = 3;
  i.get(0, l, ind);
  BOOST_CHECK_EQUAL(l, l2);
  BOOST_CHECK_EQUAL(ind, ind2);
}

BOOST_AUTO_TEST_CASE(testLevelZero) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regular(s, 2);

  HashGridIndex i(1);

  i.set(0, 1, 1);
  s.left_levelzero(&i, 0);

  HashGridIndex::level_type l, l2 = 0;
  HashGridIndex::index_type ind, ind2 = 0;

  i.get(0, l, ind);
  BOOST_CHECK_EQUAL(l, l2);
  BOOST_CHECK_EQUAL(ind, ind2);

  ind2 = 1;
  i.set(0, 1, 1);
  s.right_levelzero(&i, 0);
  i.get(0, l, ind);
  BOOST_CHECK_EQUAL(l, l2);
  BOOST_CHECK_EQUAL(ind, ind2);
}

BOOST_AUTO_TEST_CASE(testTop) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regular(s, 2);

  HashGridIndex i(1);

  i.set(0, 1, 1);
  s.left_child(&i, 0);

  s.top(&i, 0);

  HashGridIndex::level_type l, l2 = 1;
  HashGridIndex::index_type ind, ind2 = 1;
  i.get(0, l, ind);

  BOOST_CHECK_EQUAL(l, l2);
  BOOST_CHECK_EQUAL(ind, ind2);
}

BOOST_AUTO_TEST_CASE(testSeq) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regular(s, 2);

  HashGridIndex i(1);

  i.set(0, 1, 1);
  s.left_child(&i, 0);

  size_t seq = s.seq(&i);
  BOOST_CHECK(!(s.end(seq)));

  s.left_child(&i, 0);

  seq = s.seq(&i);
  BOOST_CHECK(s.end(seq));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestHashGenerator)

BOOST_AUTO_TEST_CASE(testPeriodic1D) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regularWithPeriodicBoundaries(s, 2);

  BOOST_CHECK_EQUAL(s.getSize(), 4U);
}

BOOST_AUTO_TEST_CASE(testPeriodic2D) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regularWithPeriodicBoundaries(s, 2);

  BOOST_CHECK_EQUAL(s.getSize(), 12U);

  HashGridStorage s2(2);
  g.regularWithPeriodicBoundaries(s2, 3);
  BOOST_CHECK_EQUAL(s2.getSize(), 32U);

  HashGridIndex i(2);
  i.set(0, 0, 0);
  i.set(1, 0, 0);
  BOOST_CHECK(s2.has_key(&i));

  i.set(1, 1, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(1, 2, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 2, 3);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 3, 5);
  BOOST_CHECK(!(s2.has_key(&i)));

  i.set(1, 1, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 1, 1);
  i.set(1, 1, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(1, 2, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 2, 3);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 3, 5);
  BOOST_CHECK(!(s2.has_key(&i)));

  i.set(1, 1, 1);
  BOOST_CHECK(s2.has_key(&i));
}

BOOST_AUTO_TEST_CASE(testPeriodic3D) {
  HashGridStorage s(3);
  HashGenerator g;

  g.regularWithPeriodicBoundaries(s, 2);

  BOOST_CHECK_EQUAL(s.getSize(), 32U);
}

BOOST_AUTO_TEST_CASE(testRegular1D) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regular(s, 2);

  BOOST_CHECK_EQUAL(s.getSize(), 3U);
}

BOOST_AUTO_TEST_CASE(testRegular2D) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regular(s, 2);

  BOOST_CHECK_EQUAL(s.getSize(), 5U);

  HashGridStorage s2(2);
  g.regular(s2, 3);

  BOOST_CHECK_EQUAL(s2.getSize(), 17U);

  HashGridIndex i(2);
  i.set(0, 1, 1);
  i.set(1, 1, 1);

  BOOST_CHECK(s2.has_key(&i));

  i.set(1, 2, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 2, 3);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 3, 5);
  BOOST_CHECK(!(s2.has_key(&i)));

  i.set(1, 1, 1);
  BOOST_CHECK(s2.has_key(&i));
}

BOOST_AUTO_TEST_CASE(testRegular3D) {
  HashGridStorage s(3);
  HashGenerator g;

  g.regular(s, 2);

  BOOST_CHECK_EQUAL(s.getSize(), 7U);
}

BOOST_AUTO_TEST_CASE(testRegularTruncatedBoundaries1D) {
  HashGridStorage s(1);
  HashGenerator g;

  g.regularWithBoundaries(s, 2, 1);

  BOOST_CHECK_EQUAL(s.getSize(), 5U);
}

BOOST_AUTO_TEST_CASE(testRegularTruncatedBoundaries2D) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regularWithBoundaries(s, 2, 1);

  BOOST_CHECK_EQUAL(s.getSize(), 21U);

  HashGridStorage s2(2);
  g.regularWithBoundaries(s2, 3, 1);

  BOOST_CHECK_EQUAL(s2.getSize(), 49U);

  HashGridIndex i(2);
  i.set(0, 1, 1);
  i.set(1, 1, 1);

  BOOST_CHECK(s2.has_key(&i));

  i.set(1, 2, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 2, 3);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 3, 5);
  BOOST_CHECK(!(s2.has_key(&i)));

  i.set(1, 1, 1);
  BOOST_CHECK(s2.has_key(&i));

  i.set(1, 0, 0);
  BOOST_CHECK(s2.has_key(&i));

  i.set(0, 0, 0);
  BOOST_CHECK(s2.has_key(&i));
}

BOOST_AUTO_TEST_CASE(testRegularTruncatedBoundaries3D) {
  HashGridStorage s(3);
  HashGenerator g;

  g.regularWithBoundaries(s, 2, 1);
  BOOST_CHECK_EQUAL(s.getSize(), 81U);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestHashRefinement)

BOOST_AUTO_TEST_CASE(testFreeRefine) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regular(s, 1);

  DataVector d(1);
  d[0] = 1.0;

  SurplusRefinementFunctor f(d);
  HashRefinement r;

  r.free_refine(s, f);

  BOOST_CHECK_EQUAL(s.getSize(), 5U);
}

BOOST_AUTO_TEST_CASE(testFreeRefineTruncatedBoundaries) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regularWithBoundaries(s, 1, 1);

  DataVector d(9);
  d[0] = 0.0;
  d[1] = 0.0;
  d[2] = 0.0;
  d[3] = 0.0;
  d[4] = 0.0;
  d[5] = 0.0;
  d[6] = 0.0;
  d[7] = 0.0;
  d[8] = 1.0;

  SurplusRefinementFunctor f(d);
  HashRefinementBoundaries r;

  r.free_refine(s, f);

  BOOST_CHECK_EQUAL(s.getSize(), 21U);
}

BOOST_AUTO_TEST_CASE(testFreeRefineTruncatedBoundaries2) {
  HashGridStorage s(2);
  HashGenerator g;

  g.regularWithBoundaries(s, 2, 0);

  DataVector d(17);

  for (unsigned int i = 0; i < d.getSize(); ++i) {
    d[i] = 0.0;
  }

  d[12] = 1.0;

  SurplusRefinementFunctor f(d);
  HashRefinementBoundaries r;

  r.free_refine(s, f);

  BOOST_CHECK_EQUAL(s.getSize(), 21U);
}

BOOST_AUTO_TEST_CASE(testSurplusFunctor) {
  HashGridStorage s(2);
  DataVector d(1);
  SurplusRefinementFunctor f(d);

  d[0] = -10.0;
  BOOST_CHECK_GT(f(s, 0), f.start());

  d[0] = 10.0;
  BOOST_CHECK_GT(f(s, 0), f.start());
}

BOOST_AUTO_TEST_SUITE_END()
