// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>


namespace sgpp {
namespace base {

DataVectorSP::DataVectorSP(size_t size) :
  size(size), unused(0), inc_elems(100) {
  // create new vector
  this->data = new float[size];
}

DataVectorSP::DataVectorSP(size_t size, float value) :
  DataVectorSP(size) {
  setAll(value);
}

DataVectorSP::DataVectorSP(const DataVectorSP& vec) :
  DataVectorSP(vec.size) {
  // copy data
  std::memcpy(this->data, vec.data, size * sizeof(float));
}

DataVectorSP::DataVectorSP(float* input, size_t size) :
  DataVectorSP(size) {
  // copy data
  std::memcpy(this->data, input, size * sizeof(float));
}

DataVectorSP::DataVectorSP(std::vector<float> input) :
  DataVectorSP(input.size()) {
  // copy data
  std::copy(input.begin(), input.end(), this->data);
}

DataVectorSP::DataVectorSP(std::vector<int> input) :
  DataVectorSP(input.size()) {
  // copy data
  int in = 0;

  for (std::vector<int>::iterator it = input.begin(); it < input.end();
       it++) {
    data[in] = static_cast<float>(*it);
    in++;
  }
}

void DataVectorSP::restructure(std::vector<size_t>& remainingIndex) {
  if (remainingIndex.size() > this->size) {
    throw sgpp::base::algorithm_exception("more indices than entries!");
  }

  float* newdata = new float[remainingIndex.size()];

  for (size_t i = 0; i < remainingIndex.size(); i++) {
    newdata[i] = this->data[remainingIndex[i]];
  }

  delete[] this->data;

  this->data = newdata;
  this->size = remainingIndex.size();
  this->unused = 0;
}

void DataVectorSP::resize(size_t size) {
  // don't do anyhing, if vector already has the correct size
  if (size == this->size) {
    return;
  }

  // create new vector
  float* newdata = new float[size];
  // copy entries of old vector
  std::memcpy(newdata, this->data, std::min(this->size, size)
              * sizeof(float));

  delete[] this->data;

  this->data = newdata;
  this->size = size;
  this->unused = 0;
}

void DataVectorSP::resizeZero(size_t size) {
  // don't do anyhing, if vector already has the correct size
  if (size == this->size) {
    return;
  }

  // create new vector
  float* newdata = new float[size];
  // copy entries of old vector
  std::memcpy(newdata, this->data, std::min(this->size, size)
              * sizeof(float));

  // set new elements to zero
  for (size_t i = std::min(this->size, size); i < size; i++) {
    newdata[i] = 0.0f;
  }

  delete[] this->data;

  this->data = newdata;
  this->size = size;
  this->unused = 0;
}

void DataVectorSP::addSize(size_t add) {
  // create new vector
  float* newdata = new float[(size + add)];
  // copy entries of old vector
  std::memcpy(newdata, this->data, this->size * sizeof(float));

  delete[] this->data;

  this->data = newdata;
  this->unused = add;
}

size_t DataVectorSP::append() {
  // enlarge, if necessary
  if (unused == 0) {
    addSize(this->inc_elems);
  }

  size_t x = size;
  size++;
  unused--;

  return x;
}

size_t DataVectorSP::append(float value) {
  size_t x = append();
  data[x] = value;
  return x;
}

void DataVectorSP::insert(size_t index, float value) {
  if (index > size) {
    throw sgpp::base::data_exception(
      "DataVectorSP::insert : index out of bounds");
  }

  append();
  std::memmove(data + index + 1, data + index,
               (size - 1 - index) * sizeof(float));
  data[index] = value;
}

void DataVectorSP::setAll(float value) {
  for (size_t i = 0; i < size; i++) {
    data[i] = value;
  }
}

void DataVectorSP::set(size_t i, float value) {
  data[i] = value;
}

void DataVectorSP::copyFrom(const DataVectorSP& vec) {
  // don't copy from yourself
  if (this == &vec) {
    return;
  }

  std::memcpy(this->data, vec.data, std::min(size, vec.size)
              * sizeof(float));
}

DataVectorSP& DataVectorSP::operator=(const DataVectorSP& vec) {
  if (this == &vec) {
    return *this;
  }

  if (size != vec.size) {
    // throw sgpp::base::data_exception("DataVectorSP::add : Dimensions do not match");
    delete[] data;
    size = vec.size;
    this->data = new float[size];
  }

  std::memcpy(this->data, vec.data, size * sizeof(float));
  return *this;
}

void DataVectorSP::add(const DataVectorSP& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception(
      "DataVectorSP::add : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] += vec.data[i];
  }
}

void DataVectorSP::sub(const DataVectorSP& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception(
      "DataVectorSP::sub : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] -= vec.data[i];
  }
}

void DataVectorSP::componentwise_mult(const DataVectorSP& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception(
      "DataVectorSP::componentwise_mult : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] *= vec.data[i];
  }
}

void DataVectorSP::componentwise_div(const DataVectorSP& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception(
      "DataVectorSP::componentwise_div : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] /= vec.data[i];
  }
}

float DataVectorSP::dotProduct(const DataVectorSP& vec) const {
  float sum = 0.0f;

  for (size_t i = 0; i < size; i++) {
    sum += data[i] * vec.data[i];
  }

  return sum;
}

void DataVectorSP::mult(float scalar) {
  for (size_t i = 0; i < size; i++) {
    data[i] *= scalar;
  }
}

void DataVectorSP::sqr() {
  for (size_t i = 0; i < size; i++) {
    data[i] = data[i] * data[i];
  }
}

void DataVectorSP::sqrt() {
  for (size_t i = 0; i < size; i++) {
    data[i] = std::sqrt(data[i]);
  }
}

void DataVectorSP::abs() {
  for (size_t i = 0; i < size; i++) {
    data[i] = std::fabs(data[i]);
  }
}

float DataVectorSP::sum() const {
  float result = 0.0f;

  for (size_t i = 0; i < size; i++) {
    result += data[i];
  }

  return result;
}

float DataVectorSP::maxNorm() const {
  float max = 0.0f;

  for (size_t i = 0; i < size; i++) {
    if (max < std::abs(data[i])) {
      max = std::abs(data[i]);
    }
  }

  return max;
}

float DataVectorSP::RMSNorm() const {
  float l2Norm;
  DataVectorSP temp(*this);

  temp.componentwise_mult(temp);
  l2Norm = temp.sum();
  l2Norm /= static_cast<float>(temp.getSize());
  l2Norm = std::sqrt(l2Norm);

  return l2Norm;
}

float DataVectorSP::l2Norm() const {
  float l2Norm;
  DataVectorSP temp(*this);

  temp.componentwise_mult(temp);
  l2Norm = temp.sum();
  l2Norm = std::sqrt(l2Norm);

  return l2Norm;
}

void DataVectorSP::partitionClasses(float threshold) {
  for (size_t i = 0; i < size; i++) {
    data[i] = data[i] > threshold ? 1.0f : -1.0f;
  }
}

void DataVectorSP::axpy(float a, DataVectorSP& x) {
  if (size != x.size) {
    return;
  }

  float* p_x = x.data;
  float* p_d = data;

  for (size_t i = 0; i < size; i++) {
    p_d[i] += a * p_x[i];
  }
}

void DataVectorSP::normalize() {
  normalize(0.0f);
}

void DataVectorSP::normalize(float border) {
  float min, max;
  minmax(&min, &max);

  float delta = (max - min) / (1 - 2 * border);

  for (size_t i = 0; i < size; i++) {
    data[i] = (data[i] - min) / delta + border;
  }
}

void DataVectorSP::toString(std::string& text) const {
  std::stringstream str;

  str << "[";

  for (size_t i = 0; i < size; i++) {
    if (i != 0) {
      str << ", ";
    }

    str << data[i];
  }

  str << "]";
  text = str.str();
}

std::string DataVectorSP::toString() const {
  std::string str;
  toString(str);
  return str;
}

float DataVectorSP::min() const {
  float min = std::numeric_limits<float>::infinity();

  for (size_t i = 0; i < size; i++) {
    if (min > data[i]) {
      min = data[i];
    }
  }

  return min;
}

float DataVectorSP::max() const {
  float max = -std::numeric_limits<float>::infinity();

  for (size_t i = 0; i < size; i++) {
    if (max < data[i]) {
      max = data[i];
    }
  }

  return max;
}

void DataVectorSP::minmax(float* min, float* max) const {
  float min_t = std::numeric_limits<float>::infinity();
  float max_t = -std::numeric_limits<float>::infinity();

  for (size_t i = 0; i < size; i++) {
    if (min_t > data[i]) {
      min_t = data[i];
    }

    if (max_t < data[i]) {
      max_t = data[i];
    }
  }

  (*min) = min_t;
  (*max) = max_t;
}

float* DataVectorSP::getPointer() {
  return data;
}

const float* DataVectorSP::getPointer() const {
  return data;
}

DataVectorSP::~DataVectorSP() {
  delete[] data;
}

size_t DataVectorSP::getNumberNonZero() const {
  size_t nonZero = 0;

  for (size_t i = 0; i < size; i++) {
    if (std::abs(data[i]) > 0.0f) {
      nonZero++;
    }
  }

  return nonZero;
}

}  // namespace base
}  // namespace sgpp
