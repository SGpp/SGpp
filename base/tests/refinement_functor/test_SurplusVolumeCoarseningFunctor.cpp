// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>

#include <cmath>

using sgpp::base::DataVector;
using sgpp::base::HashGenerator;
using sgpp::base::HashCoarsening;
using sgpp::base::HashGridStorage;
using sgpp::base::SurplusVolumeCoarseningFunctor;

BOOST_AUTO_TEST_SUITE(TestSurplusVolumeCoarseningFunctor)

BOOST_AUTO_TEST_CASE(testCoarsen) {
  HashGridStorage storage(2);
  HashGenerator generator;
  HashCoarsening coarsen;

  generator.regular(storage, 2);

  DataVector alpha(storage.getSize(), 0.5);
  alpha[0] = 1.0;  // these surpluses corresponds to constant function 1

  // this should remove all children, i.e. the whole level 2
  SurplusVolumeCoarseningFunctor functor(alpha, storage.getSize(), 0.0625);
  coarsen.free_coarsen(storage, functor, alpha);

  BOOST_CHECK_EQUAL(storage.getSize(), 1);
}

BOOST_AUTO_TEST_SUITE_END()
