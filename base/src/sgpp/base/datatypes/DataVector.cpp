// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

DataVector::DataVector() : DataVector(0) {}

DataVector::DataVector(size_t size) : DataVector(size, 0.0) {}

DataVector::DataVector(size_t size, double value) { this->assign(size, value); }

DataVector::DataVector(double* input, size_t size)
    : std::vector<double>(input, input + size) {}

DataVector::DataVector(std::vector<double> input)
    : std::vector<double>(input) {}

DataVector::DataVector(std::initializer_list<double> input)
    : std::vector<double>(input) {}

DataVector::DataVector(std::vector<int> input) {
  // copy data
  this->reserve(input.size());
  for (auto iter = input.begin(); iter < input.end(); ++iter) {
    this->emplace_back(static_cast<double>(*iter));
  }
}

DataVector DataVector::fromFile(const std::string& fileName) {
  std::ifstream f(fileName, std::ifstream::in);
  f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  std::string content;
  content.assign(std::istreambuf_iterator<char>(f),
                 std::istreambuf_iterator<char>());
  return DataVector::fromString(content);
}

DataVector DataVector::fromString(const std::string& serializedVector) {
  DataVector v;

  enum class PARSER_STATE { INIT, VALUE, COMMAEND, END };

  PARSER_STATE state = PARSER_STATE::INIT;

  size_t i = 0;
  while (i < serializedVector.size()) {
    char c = serializedVector[i];

    if (std::isspace(c)) {
      i += 1;
      continue;
    }

    if (state == PARSER_STATE::INIT) {
      if (c != '[') {
        throw;
      }
      state = PARSER_STATE::VALUE;
      i += 1;
    } else if (state == PARSER_STATE::VALUE) {
      //      size_t next;
      //      double value = std::atof(&(serializedVector[i]));
      size_t endNumber = i;
      while (serializedVector[endNumber] != ',' &&
             serializedVector[endNumber] != ']')
        ++endNumber;
      std::stringstream stream;
      for (size_t j = i; j < endNumber; ++j) {
        stream << serializedVector[j];
      }
      std::string shortString(stream.str());
      //      double value = std::atof(shortString.c_str());
      double value = std::stod(shortString);
      v.append(value);
      state = PARSER_STATE::COMMAEND;
      //      i += next;
      //      while (serializedVector[i] != ',' && serializedVector[i] != ']')
      //      ++i;
      i = endNumber;

    } else if (state == PARSER_STATE::COMMAEND) {
      if (c == ',') {
        state = PARSER_STATE::VALUE;
        ++i;
      } else if (c == ']') {
        state = PARSER_STATE::END;
        ++i;
      }
    } else if (state == PARSER_STATE::END) {
      // only reached if a non-whitespace character was encountered after
      // closing brace
      throw data_exception("error: could not parse DataVector file");
    }
  }
  return v;
}

void DataVector::resizeZero(size_t size) { this->resize(size); }

void DataVector::restructure(std::vector<size_t>& remainingIndex) {
  DataVector oldVector(*this);
  this->resize(remainingIndex.size());

  for (size_t i = 0; i < remainingIndex.size(); ++i) {
    (*this)[i] = oldVector[remainingIndex[i]];
  }
}

void DataVector::remove(std::vector<size_t>& indexesToRemove) {
  DataVector oldVector(*this);
  std::vector<bool> willBeRemoved(this->size(), false);

  // Count the indexes to remove for the case when there are duplicates in
  // indexesToRemove
  size_t numIndexesToRemove = 0;
  for (size_t idx : indexesToRemove) {
    if (!willBeRemoved[idx]) {
      willBeRemoved[idx] = true;
      ++numIndexesToRemove;
    }
  }

  this->resize(oldVector.size() - numIndexesToRemove);
  size_t j = 0;
  for (size_t i = 0; i < oldVector.size(); ++i) {
    if (!willBeRemoved[i]) {
      (*this)[j] = oldVector[i];
      ++j;
    }
  }
}

size_t DataVector::append() { return this->append(0.0); }

size_t DataVector::append(double value) {
  this->emplace_back(value);
  return this->size() - 1;
}

void DataVector::append(DataVector::iterator first, DataVector::iterator last) {
  this->insert(this->end(), first, last);
}

void DataVector::insert(size_t index, double value) {
  if (index > this->size()) {
    throw sgpp::base::data_exception(
        "DataVector::insert : index out of bounds");
  }

  this->insert(this->begin() + index, value);
}

void DataVector::setAll(double value) {
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = value;
  }
}

void DataVector::set(size_t i, double value) { (*this)[i] = value; }

void DataVector::copyFrom(const DataVector& vec) {
  if (*this == vec) {
    return;
  }
  std::copy(vec.begin(), vec.begin() + std::min(this->size(), vec.size()),
            this->begin());
}

void DataVector::add(const DataVector& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception(
        "DataVector::add : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] += vec[i];
  }
}

void DataVector::accumulate(const DataVector& vec) {
  if (this->correction.size() != this->size()) {
    this->correction.resize(this->size());
  }

  double y, t;

  for (size_t i = 0; i < this->size(); ++i) {
    y = vec[i] - this->correction[i];
    t = (*this)[i] + y;
    this->correction[i] = (t - (*this)[i]) - y;
    (*this)[i] = t;
  }
}

void DataVector::sub(const DataVector& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception(
        "DataVector::sub : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] -= vec[i];
  }
}

void DataVector::componentwise_mult(const DataVector& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception(
        "DataVector::componentwise_mult : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] *= vec[i];
  }
}

void DataVector::componentwise_div(const DataVector& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception(
        "DataVector::componentwise_div : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] /= vec[i];
  }
}

double DataVector::dotProduct(const DataVector& vec) const {
  double sum = 0.0;

  for (size_t i = 0; i < this->size(); ++i) {
    sum += (*this)[i] * vec[i];
  }

  return sum;
}

void DataVector::mult(double scalar) {
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] *= scalar;
  }
}

void DataVector::sqr() {
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = (*this)[i] * (*this)[i];
  }
}

void DataVector::sqrt() {
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = std::sqrt((*this)[i]);
  }
}

void DataVector::abs() {
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = std::abs((*this)[i]);
  }
}

double DataVector::sum() const {
  double result = 0.0;

  for (size_t i = 0; i < this->size(); ++i) {
    result += (*this)[i];
  }

  return result;
}

double DataVector::maxNorm() const {
  double max = 0.0;

  for (size_t i = 0; i < this->size(); ++i) {
    if (max < std::abs((*this)[i])) {
      max = std::abs((*this)[i]);
    }
  }

  return max;
}

double DataVector::RMSNorm() const {
  double rmsNorm;
  DataVector temp(*this);

  temp.sqr();
  rmsNorm = temp.sum();
  rmsNorm /= static_cast<double>(temp.getSize());
  rmsNorm = std::sqrt(rmsNorm);

  return rmsNorm;
}

double DataVector::l2Norm() const {
  double l2Norm;
  DataVector temp(*this);

  temp.componentwise_mult(temp);
  l2Norm = temp.sum();
  l2Norm = std::sqrt(l2Norm);

  return l2Norm;
}

double DataVector::min() const {
  double min = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < this->size(); ++i) {
    if (min > (*this)[i]) {
      min = (*this)[i];
    }
  }

  return min;
}

double DataVector::max() const {
  double max = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < this->size(); ++i) {
    if (max < (*this)[i]) {
      max = (*this)[i];
    }
  }

  return max;
}

void DataVector::minmax(double* min, double* max) const {
  double min_t = std::numeric_limits<double>::infinity();
  double max_t = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < this->size(); ++i) {
    if (min_t > (*this)[i]) {
      min_t = (*this)[i];
    }

    if (max_t < (*this)[i]) {
      max_t = (*this)[i];
    }
  }

  (*min) = min_t;
  (*max) = max_t;
}

void DataVector::axpy(double a, DataVector& x) {
  if (this->size() != x.size()) {
    return;
  }

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] += a * x[i];
  }
}

double* DataVector::getPointer() { return const_cast<double*>(this->data()); }

const double* DataVector::getPointer() const { return this->data(); }

size_t DataVector::getNumberNonZero() const {
  size_t nonZero = 0;

  for (size_t i = 0; i < this->size(); ++i) {
    if (std::abs((*this)[i]) > 0.0) {
      ++nonZero;
    }
  }

  return nonZero;
}

void DataVector::partitionClasses(double threshold) {
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = (*this)[i] > threshold ? 1.0 : -1.0;
  }
}

void DataVector::normalize() { normalize(0.0); }

void DataVector::normalize(double border) {
  double min, max;
  minmax(&min, &max);

  double delta = (max - min) / (1.0 - 2.0 * border);

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = ((*this)[i] - min) / delta + border;
  }
}

void DataVector::toString(std::string& text) const {
  std::stringstream str;

  str << std::scientific;
  str.precision(20);

  str << "[";

  for (size_t i = 0; i < this->size(); ++i) {
    if (i != 0) {
      str << ", ";
      // add linebreaks for readability
      if (i % 20 == 0) {
        str << std::endl;
      }
    }

    str << (*this)[i];
  }

  str << "]";
  text = str.str();
}

std::string DataVector::toString() const {
  std::string str;
  toString(str);
  return str;
}

void DataVector::toFile(const std::string& fileName) const {
  std::ofstream f(fileName, std::ofstream::out);
  f << this->toString();
  f.close();
}

}  // namespace base
}  // namespace sgpp
