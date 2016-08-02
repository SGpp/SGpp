/*
 * FunctionLookupTable.hpp
 *
 *  Created on: 07.11.2015
 *      Author: david
 */

#ifndef FUNCTIONLOOKUPTABLE_HPP_
#define FUNCTIONLOOKUPTABLE_HPP_

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/MultiFunction.hpp>

#include <vector>
#include <memory>
#include <unordered_map>
#include <string>
#include <mutex>

namespace SGPP {
namespace combigrid {

class DataVectorEqualTo {
public:
	bool operator()(base::DataVector const &first, base::DataVector const &second) const {
		if (first.getSize() != second.getSize()) {
			return false;
		}

		for (size_t i = 0; i < first.getSize(); ++i) {
			if (first[i] != second[i]) {
				return false;
			}
		}

		return true;
	}
};

class DataVectorHash {
public:
	size_t operator()(base::DataVector const &vec) const {
		std::hash<SGPP::float_t> h;
		size_t result = 0;
		for (size_t i = 0; i < vec.getSize(); ++i) {
			result ^= h(vec[i]) * (i + 1);
		}
		return result;
	}
};

class FunctionLookupTable {
	std::shared_ptr<std::unordered_map<base::DataVector, SGPP::float_t, DataVectorHash, DataVectorEqualTo>> hashmap;
	MultiFunction func;
	std::mutex tableMutex;

public:
	FunctionLookupTable(MultiFunction func);

	SGPP::float_t operator()(base::DataVector const &x);

	SGPP::float_t eval(base::DataVector const &x);
	SGPP::float_t evalThreadsafe(base::DataVector const &x);

	bool containsEntry(base::DataVector const &x);

	void addEntry(base::DataVector const &x, SGPP::float_t y);

	std::string serialize();

	void deserialize(std::string const &value);

	size_t getNumEntries() const;
};

}
} /* namespace utils */

#endif /* FUNCTIONLOOKUPTABLE_HPP_ */
