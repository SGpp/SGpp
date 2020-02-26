// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridPoint.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp>

using sgpp::base::DataVector;
using sgpp::base::HashGenerator;
using sgpp::base::HashGridPoint;
using sgpp::base::HashGridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::SubspaceRefinement;
using sgpp::base::SurplusRefinementFunctor;

BOOST_AUTO_TEST_SUITE(TestSubspaceRefinement)

BOOST_AUTO_TEST_CASE(testFreeRefineTrivial) {
  HashGridStorage storage(2);
  HashGenerator generator;

  generator.regular(storage, 1);

  DataVector datavector(1);
  datavector[0] = 1.0;

  SurplusRefinementFunctor functor(datavector);
  HashRefinement* hash_refinement = new HashRefinement();
  SubspaceRefinement subspace_refinement(hash_refinement);

  subspace_refinement.free_refine(storage, functor);

  BOOST_CHECK_EQUAL(storage.getSize(), 5);

  delete hash_refinement;
}

BOOST_AUTO_TEST_CASE(testFreeRefineSubspaceAnisotropic) {
  HashGridStorage storage(2);
  HashGenerator generator;

  generator.regular(storage, 3);

  DataVector data_vector(storage.getSize());
  data_vector[9] = 2.0;

  SurplusRefinementFunctor functor(data_vector, 1);
  HashRefinement* hash_refinement = new HashRefinement();

  SubspaceRefinement subspace_refinement(hash_refinement);

  subspace_refinement.free_refine(storage, functor);

  BOOST_CHECK_EQUAL(storage.getSize(), 33);

  for (size_t i = 0; i < storage.getSize(); i++) {
    HashGridPoint& index = storage.getPoint(i);
    BOOST_CHECK((index.getIndex(0) == 4) == false);
  }

  delete hash_refinement;
}

BOOST_AUTO_TEST_CASE(testFreeRefineSubspaceIsotropic) {
  HashGridStorage storage(2);
  HashGenerator generator;

  generator.regular(storage, 3);

  DataVector data_vector(storage.getSize());
  data_vector[13] = 2.0;

  SurplusRefinementFunctor functor(data_vector, 1);
  HashRefinement* hash_refinement = new HashRefinement();

  SubspaceRefinement subspace_refinement(hash_refinement);

  subspace_refinement.free_refine(storage, functor);

  BOOST_CHECK_EQUAL(storage.getSize(), 33);

  for (size_t i = 0; i < storage.getSize(); i++) {
    HashGridPoint& index = storage.getPoint(i);
    BOOST_CHECK((index.getIndex(0) == 4 || index.getIndex(1) == 4) == false);
  }

  delete hash_refinement;
}

BOOST_AUTO_TEST_SUITE_END()
