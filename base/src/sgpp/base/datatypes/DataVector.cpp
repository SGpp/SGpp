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
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace sgpp {
namespace base {

DataVector::DataVector() : correction(NULL), size(0), unused(0), inc_elems(100) {
  // create new vector
  this->data = new double[0];
}

DataVector::DataVector(size_t size) : correction(NULL), size(size), unused(0), inc_elems(100) {
  // create new vector
  this->data = new double[size];
}

DataVector::DataVector(size_t size, double value) : DataVector(size) { setAll(value); }

DataVector::DataVector(const DataVector& vec) : DataVector(vec.size) {
  // copy data
  std::memcpy(this->data, vec.data, size * sizeof(double));
}

DataVector::DataVector(double* input, size_t size) : DataVector(size) {
  // copy data
  std::memcpy(this->data, input, size * sizeof(double));
}

DataVector::DataVector(std::vector<double> input) : DataVector(input.size()) {
  // copy data
  std::copy(input.begin(), input.end(), this->data);
}

DataVector::DataVector(std::vector<int> input) : DataVector(input.size()) {
  // copy data
  int in = 0;

  for (std::vector<int>::iterator it = input.begin(); it < input.end(); it++) {
    data[in] = static_cast<double>(*it);
    in++;
  }
}

DataVector::DataVector(DataVectorDefinition& DataVectorDef) {
  setDataVectorDefinition(DataVectorDef);
}

DataVector DataVector::fromFile(const std::string& fileName) {
  std::ifstream f(fileName, std::ifstream::in);
  f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  std::string content;
  content.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
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
      double value = std::atof(&(serializedVector[i]));
      v.append(value);
      state = PARSER_STATE::COMMAEND;
      //      i += next;
      while (serializedVector[i] != ',' && serializedVector[i] != ']') i++;

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
      throw data_exception("error: could not parse DataVector file");
    }
  }
  return v;
}

void DataVector::getDataVectorDefinition(DataVectorDefinition& DataVectorDef) {
  DataVectorDef.size = this->size;
  DataVectorDef.unused = this->unused;
  DataVectorDef.inc_elems = this->inc_elems;
  DataVectorDef.data = this->data;
}

void DataVector::setDataVectorDefinition(DataVectorDefinition& DataVectorDef) {
  this->size = DataVectorDef.size;
  this->unused = DataVectorDef.unused;
  this->inc_elems = DataVectorDef.inc_elems;
  this->data = DataVectorDef.data;
}

void DataVector::restructure(std::vector<size_t>& remainingIndex) {
  if (remainingIndex.size() > this->size) {
    throw sgpp::base::algorithm_exception("more indices than entries!");
  }

  double* newdata = new double[remainingIndex.size()];

  for (size_t i = 0; i < remainingIndex.size(); i++) {
    newdata[i] = this->data[remainingIndex[i]];
  }

  delete[] this->data;

  this->data = newdata;
  this->size = remainingIndex.size();
  this->unused = 0;
}

void DataVector::resize(size_t size) {
  // don't do anyhing, if vector already has the correct size
  if (size == this->size) {
    return;
  }

  // create new vector
  double* newdata = new double[size];
  // copy entries of old vector
  std::memcpy(newdata, this->data, std::min(this->size, size) * sizeof(double));
  delete[] this->data;
  this->data = newdata;

  if (this->correction != NULL) {
    double* newcorrection = new double[size];
    std::memcpy(newcorrection, this->correction, std::min(this->size, size) * sizeof(double));
    delete[] this->correction;
    this->correction = newcorrection;
  }

  this->size = size;
  this->unused = 0;
}

void DataVector::resizeZero(size_t size) {
  // don't do anyhing, if vector already has the correct size
  if (size == this->size) {
    return;
  }

  // create new vector
  double* newdata = new double[size];
  // copy entries of old vector
  std::memcpy(newdata, this->data, std::min(this->size, size) * sizeof(double));

  // set new elements to zero
  for (size_t i = std::min(this->size, size); i < size; i++) {
    newdata[i] = 0.0;
  }

  delete[] this->data;

  this->data = newdata;
  this->size = size;
  this->unused = 0;
}

void DataVector::addSize(size_t add) {
  // create new vector
  double* newdata = new double[(size + add)];
  // copy entries of old vector
  std::memcpy(newdata, this->data, this->size * sizeof(double));

  delete[] this->data;

  this->data = newdata;
  this->unused = add;
}

size_t DataVector::append() {
  // enlarge, if necessary
  if (unused == 0) {
    addSize(this->inc_elems);
  }

  size_t x = size;
  size++;
  unused--;

  return x;
}

size_t DataVector::append(double value) {
  size_t x = append();
  data[x] = value;
  return x;
}

void DataVector::insert(size_t index, double value) {
  if (index > size) {
    throw sgpp::base::data_exception("DataVector::insert : index out of bounds");
  }

  append();
  std::memmove(data + index + 1, data + index, (size - 1 - index) * sizeof(double));
  data[index] = value;
}

/*
 int DataVector::addValue() {
 if (unused == 0) {
 addSize(size);
 }

 int x = size;

 size++;
 unused--;

 return x;
 }
 */

void DataVector::setAll(double value) {
  for (size_t i = 0; i < size; i++) {
    data[i] = value;
  }
}

void DataVector::set(size_t i, double value) { data[i] = value; }

void DataVector::copyFrom(const DataVector& vec) {
  // don't copy from yourself
  if (this == &vec) {
    return;
  }

  /*
   if (size != vec.size) {
   delete[] data;
   size = vec.size;
   this->data = new double[size];
   }
   */
  std::memcpy(this->data, vec.data, std::min(size, vec.size) * sizeof(double));
}

DataVector& DataVector::operator=(const DataVector& vec) {
  if (this == &vec) {
    return *this;
  }

  if (size != vec.size) {
    throw sgpp::base::data_exception("DataVector::add : Dimensions do not match");
    //        delete[] data;
    //        size = vec.size;
    //        this->data = new double[size];
  }

  std::memcpy(this->data, vec.data, size * sizeof(double));
  return *this;
}

void DataVector::add(const DataVector& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception("DataVector::add : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] += vec.data[i];
  }
}

void DataVector::accumulate(const DataVector& vec) {
  if (this->correction == NULL) {
    this->correction = new double[size];
    std::memset(this->correction, 0, size * sizeof(double));
  }

  double y, t;

  for (size_t i = 0; i < size; i++) {
    y = vec[i] - this->correction[i];
    t = data[i] + y;
    this->correction[i] = (t - data[i]) - y;
    data[i] = t;
  }
}

void DataVector::sub(const DataVector& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception("DataVector::sub : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] -= vec.data[i];
  }
}

void DataVector::componentwise_mult(const DataVector& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception("DataVector::componentwise_mult : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] *= vec.data[i];
  }
}

void DataVector::componentwise_div(const DataVector& vec) {
  if (size != vec.size) {
    throw sgpp::base::data_exception("DataVector::componentwise_div : Dimensions do not match");
  }

  for (size_t i = 0; i < size; i++) {
    data[i] /= vec.data[i];
  }
}

double DataVector::dotProduct(const DataVector& vec) const {
  double sum = 0.0;

  for (size_t i = 0; i < size; i++) {
    sum += data[i] * vec.data[i];
  }

  return sum;
}

void DataVector::mult(double scalar) {
  for (size_t i = 0; i < size; i++) {
    data[i] *= scalar;
  }
}

void DataVector::sqr() {
  for (size_t i = 0; i < size; i++) {
    data[i] = data[i] * data[i];
  }
}

void DataVector::sqrt() {
  for (size_t i = 0; i < size; i++) {
    data[i] = std::sqrt(data[i]);
  }
}

void DataVector::abs() {
  for (size_t i = 0; i < size; i++) {
    data[i] = std::fabs(data[i]);
  }
}

double DataVector::sum() const {
  double result = 0.0;

  for (size_t i = 0; i < size; i++) {
    result += data[i];
  }

  return result;
}

double DataVector::maxNorm() const {
  double max = 0.0;

  for (size_t i = 0; i < size; i++) {
    if (max < fabs(data[i])) {
      max = fabs(data[i]);
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

void DataVector::partitionClasses(double threshold) {
  for (size_t i = 0; i < size; i++) {
    data[i] = data[i] > threshold ? 1.0 : -1.0;
  }
}

void DataVector::axpy(double a, DataVector& x) {
  if (size != x.size) {
    return;
  }

  double* p_x = x.data;
  double* p_d = data;

  for (size_t i = 0; i < size; i++) {
    p_d[i] += a * p_x[i];
  }
}

void DataVector::normalize() { normalize(0.0); }

void DataVector::normalize(double border) {
  double min, max;
  minmax(&min, &max);

  double delta = (max - min) / (1 - 2 * border);

  for (size_t i = 0; i < size; i++) {
    data[i] = (data[i] - min) / delta + border;
  }
}

void DataVector::toString(std::string& text) const {
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

double DataVector::min() const {
  double min = INFINITY;

  for (size_t i = 0; i < size; i++) {
    if (min > data[i]) {
      min = data[i];
    }
  }

  return min;
}

double DataVector::max() const {
  double max = -INFINITY;

  for (size_t i = 0; i < size; i++) {
    if (max < data[i]) {
      max = data[i];
    }
  }

  return max;
}

void DataVector::minmax(double* min, double* max) const {
  double min_t = INFINITY;
  double max_t = -INFINITY;

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

double* DataVector::getPointer() { return data; }

const double* DataVector::getPointer() const { return data; }

DataVector::~DataVector() {
  if (this->correction != NULL) {
    delete[] correction;
  }

  delete[] data;
}

size_t DataVector::getNumberNonZero() const {
  size_t nonZero = 0;

  for (size_t i = 0; i < size; i++) {
    if (fabs(data[i]) > 0.0) {
      nonZero++;
    }
  }

  return nonZero;
}

}  // namespace base
}  // namespace sgpp
