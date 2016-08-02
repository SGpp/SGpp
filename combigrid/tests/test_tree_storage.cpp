/*
 * test_tree_storage.cpp
 *
 *  Created on: 16.12.2015
 *      Author: david
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <iostream>

using namespace SGPP::combigrid;

BOOST_AUTO_TEST_CASE(testTreeStorageGetSet) {
	TreeStorage<int> storage(3);

	MultiIndex index(3, 0);

	BOOST_CHECK(!storage.containsIndex(index));

	BOOST_CHECK_EQUAL(storage.get(index), 0);

	BOOST_CHECK(storage.containsIndex(index));

	storage.get(index) = 5;

	BOOST_CHECK_EQUAL(storage.get(index), 5);

	index[1] = 5;
	index[2] = 5;

	storage.set(index, 7);

	BOOST_CHECK_EQUAL(storage.get(index), 7);
}

BOOST_AUTO_TEST_CASE(testTreeStorageDataIterator) {
	TreeStorage<int> storage(3);

	MultiIndex index(3, 0);
	storage.set(index, 1);

	index[2] = 2;
	storage.set(index, 2);

	index[0] = 2;
	storage.set(index, 3);

	auto it = storage.getStoredDataIterator();

	BOOST_CHECK(it->isValid());

	BOOST_CHECK_EQUAL(it->indexAt(0), 0);
	BOOST_CHECK_EQUAL(it->indexAt(1), 0);
	BOOST_CHECK_EQUAL(it->indexAt(2), 0);

	BOOST_CHECK_EQUAL(it->value(), 1);

	BOOST_CHECK_EQUAL(it->moveToNext(), 0);
	BOOST_CHECK_EQUAL(it->value(), 2);

	//BOOST_CHECK_EQUAL(it->moveToNext(), 2);
	//BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 2);
	BOOST_CHECK_EQUAL(it->value(), 3);

	BOOST_CHECK_EQUAL(it->moveToNext(), -1);

	BOOST_CHECK(!it->isValid());
}

BOOST_AUTO_TEST_CASE(testTreeStorageGuidedIterator) {
	TreeStorage<int> storage(3);

	MultiIndex bounds(3, 2);
	MultiIndexIterator multiIter(bounds);
	MultiIndex index(3, 0);
	storage.set(index, 1);

	index[2] = 2;
	storage.set(index, 2);

	index[0] = 2;
	storage.set(index, 3);

	auto it = storage.getGuidedIterator(multiIter);

	BOOST_CHECK(it->isValid());

	BOOST_CHECK_EQUAL(it->indexAt(0), 0);
	BOOST_CHECK_EQUAL(it->indexAt(1), 0);
	BOOST_CHECK_EQUAL(it->indexAt(2), 0);

	BOOST_CHECK_EQUAL(it->value(), 1);

	BOOST_CHECK_EQUAL(it->moveToNext(), 0);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 1);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 0);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 2);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 0);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 1);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), 0);
	BOOST_CHECK_EQUAL(it->value(), 0);

	BOOST_CHECK_EQUAL(it->moveToNext(), -1);
	BOOST_CHECK_EQUAL(it->isValid(), false);
}
