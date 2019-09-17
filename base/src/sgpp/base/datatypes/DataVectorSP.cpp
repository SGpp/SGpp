// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace sgpp {
namespace base {

DataVectorSP::DataVectorSP() : DataVectorSP(0) {}

DataVectorSP::DataVectorSP(size_t size) : DataVectorSP(size, 0.0) {}

DataVectorSP::DataVectorSP(size_t size, float value) { this->assign(size, value); }

DataVectorSP::DataVectorSP(float* input, size_t size) : std::vector<float>(input, input + size) {}

DataVectorSP::DataVectorSP(std::vector<float> input) : std::vector<float>(input) {}

DataVectorSP::DataVectorSP(std::initializer_list<float> input) : std::vector<float>(input) {}

DataVectorSP::DataVectorSP(std::vector<int> input) {
  // copy data
  this->reserve(input.size());
  for (auto iter = input.begin(); iter < input.end(); ++iter) {
    this->emplace_back(static_cast<float>(*iter));
  }
}

DataVectorSP DataVectorSP::fromFile(const std::string& fileName) {
  std::ifstream f(fileName, std::ifstream::in);
  f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  std::string content;
  content.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
  return DataVectorSP::fromString(content);
}

DataVectorSP DataVectorSP::fromString(const std::string& serializedVector) {
  DataVectorSP v;

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
      //      float value = std::atof(&(serializedVector[i]));
      size_t endNumber = i;
      while (serializedVector[endNumber] != ',' && serializedVector[endNumber] != ']') endNumber++;
      std::stringstream stream;
      for (size_t j = i; j < endNumber; j++) {
        stream << serializedVector[j];
      }
      std::string shortString(stream.str());
      //      float value = std::atof(shortString.c_str());
      float value = std::stof(shortString);
      v.append(value);
      state = PARSER_STATE::COMMAEND;
      //      i += next;
      //      while (serializedVector[i] != ',' && serializedVector[i] != ']') i++;
      i = endNumber;

    } else if (state == PARSER_STATE::COMMAEND) {
      if (c == ',') {
        state = PARSER_STATE::VALUE;
        i++;
      } else if (c == ']') {
        state = PARSER_STATE::END;
        i++;
      }
    } else if (state == PARSER_STATE::END) {
      // only reached if a non-whitespace character was encountered after closing brace
      throw data_exception("error: could not parse DataVectorSP file");
    }
  }
  return v;
}

void DataVectorSP::resizeZero(size_t size) { this->resize(size); }

void DataVectorSP::restructure(std::vector<size_t>& remainingIndex) {
  DataVectorSP oldVector(*this);
  this->resize(remainingIndex.size());

  for (size_t i = 0; i < remainingIndex.size(); i++) {
    (*this)[i] = oldVector[remainingIndex[i]];
  }
}

void DataVectorSP::remove(std::vector<size_t>& indexesToRemove) {
  DataVectorSP oldVector(*this);
  std::vector<bool> willBeRemoved(this->size(), false);

  // Count the indexes to remove for the case when there are duplicates in indexesToRemove
  size_t numIndexesToRemove = 0;
  for (size_t i = 0; i < oldVector.size(); i++) {
    size_t idx = indexesToRemove[i];
    if (!willBeRemoved[idx]) {
      willBeRemoved[idx] = true;
      numIndexesToRemove++;
    }
  }

  this->resize(oldVector.size() - numIndexesToRemove);
  size_t j = 0;
  for (size_t i = 0; i < oldVector.size(); i++) {
    if (!willBeRemoved[i]) {
      (*this)[j] = oldVector[i];
      j++;
    }
  }
}

size_t DataVectorSP::append() { return this->append(0.0f); }

size_t DataVectorSP::append(float value) {
  this->emplace_back(value);
  return this->size() - 1;
}

void DataVectorSP::insert(size_t index, float value) {
  if (index > this->size()) {
    throw sgpp::base::data_exception("DataVectorSP::insert : index out of bounds");
  }

  this->insert(this->begin() + index, value);
}

void DataVectorSP::setAll(float value) {
  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] = value;
  }
}

void DataVectorSP::set(size_t i, float value) { (*this)[i] = value; }

void DataVectorSP::copyFrom(const DataVectorSP& vec) {
  if (*this == vec) {
    return;
  }
  std::copy(vec.begin(), vec.begin() + std::min(this->size(), vec.size()), this->begin());
}

void DataVectorSP::add(const DataVectorSP& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception("DataVectorSP::add : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] += vec[i];
  }
}

void DataVectorSP::accumulate(const DataVectorSP& vec) {
  if (this->correction.size() != this->size()) {
    this->correction.resize(this->size());
  }

  float y, t;

  for (size_t i = 0; i < this->size(); i++) {
    y = vec[i] - this->correction[i];
    t = (*this)[i] + y;
    this->correction[i] = (t - (*this)[i]) - y;
    (*this)[i] = t;
  }
}

void DataVectorSP::sub(const DataVectorSP& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception("DataVectorSP::sub : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] -= vec[i];
  }
}

void DataVectorSP::componentwise_mult(const DataVectorSP& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception("DataVectorSP::componentwise_mult : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] *= vec[i];
  }
}

void DataVectorSP::componentwise_div(const DataVectorSP& vec) {
  if (this->size() != vec.size()) {
    throw sgpp::base::data_exception("DataVectorSP::componentwise_div : Dimensions do not match");
  }

  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] /= vec[i];
  }
}

float DataVectorSP::dotProduct(const DataVectorSP& vec) const {
  float sum = 0.0f;

  for (size_t i = 0; i < this->size(); i++) {
    sum += (*this)[i] * vec[i];
  }

  return sum;
}

void DataVectorSP::mult(float scalar) {
  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] *= scalar;
  }
}

void DataVectorSP::sqr() {
  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] = (*this)[i] * (*this)[i];
  }
}

void DataVectorSP::sqrt() {
  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] = std::sqrt((*this)[i]);
  }
}

void DataVectorSP::abs() {
  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] = std::abs((*this)[i]);
  }
}

float DataVectorSP::sum() const {
  float result = 0.0f;

  for (size_t i = 0; i < this->size(); i++) {
    result += (*this)[i];
  }

  return result;
}

float DataVectorSP::maxNorm() const {
  float max = 0.0f;

  for (size_t i = 0; i < this->size(); i++) {
    if (max < std::abs((*this)[i])) {
      max = std::abs((*this)[i]);
    }
  }

  return max;
}

float DataVectorSP::RMSNorm() const {
  float rmsNorm;
  DataVectorSP temp(*this);

  temp.sqr();
  rmsNorm = temp.sum();
  rmsNorm /= static_cast<float>(temp.getSize());
  rmsNorm = std::sqrt(rmsNorm);

  return rmsNorm;
}

float DataVectorSP::l2Norm() const {
  float l2Norm;
  DataVectorSP temp(*this);

  temp.componentwise_mult(temp);
  l2Norm = temp.sum();
  l2Norm = std::sqrt(l2Norm);

  return l2Norm;
}

float DataVectorSP::min() const {
  float min = std::numeric_limits<float>::infinity();

  for (size_t i = 0; i < this->size(); i++) {
    if (min > (*this)[i]) {
      min = (*this)[i];
    }
  }

  return min;
}

float DataVectorSP::max() const {
  float max = -std::numeric_limits<float>::infinity();

  for (size_t i = 0; i < this->size(); i++) {
    if (max < (*this)[i]) {
      max = (*this)[i];
    }
  }

  return max;
}

void DataVectorSP::minmax(float* min, float* max) const {
  float min_t = std::numeric_limits<float>::infinity();
  float max_t = -std::numeric_limits<float>::infinity();

  for (size_t i = 0; i < this->size(); i++) {
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

void DataVectorSP::axpy(float a, DataVectorSP& x) {
  if (this->size() != x.size()) {
    return;
  }

  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] += a * x[i];
  }
}

float* DataVectorSP::getPointer() { return const_cast<float*>(this->data()); }

const float* DataVectorSP::getPointer() const { return this->data(); }

size_t DataVectorSP::getNumberNonZero() const {
  size_t nonZero = 0;

  for (size_t i = 0; i < this->size(); i++) {
    if (std::abs((*this)[i]) > 0.0f) {
      nonZero++;
    }
  }

  return nonZero;
}

void DataVectorSP::partitionClasses(float threshold) {
  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] = (*this)[i] > threshold ? 1.0f : -1.0f;
  }
}

void DataVectorSP::normalize() { normalize(0.0f); }

void DataVectorSP::normalize(float border) {
  float min, max;
  minmax(&min, &max);

  float delta = (max - min) / (1.0f - 2.0f * border);

  for (size_t i = 0; i < this->size(); i++) {
    (*this)[i] = ((*this)[i] - min) / delta + border;
  }
}

void DataVectorSP::toString(std::string& text) const {
  std::stringstream str;

  str << std::scientific;
  str.precision(20);

  str << "[";

  for (size_t i = 0; i < this->size(); i++) {
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

std::string DataVectorSP::toString() const {
  std::string str;
  toString(str);
  return str;
}

void DataVectorSP::toFile(const std::string& fileName) const {
  std::ofstream f(fileName, std::ofstream::out);
  f << this->toString();
  f.close();
}

}  // namespace base
}  // namespace sgpp
