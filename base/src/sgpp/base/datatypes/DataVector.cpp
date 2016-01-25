// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>
//#include <sgpp/base/tools/AlignedMemory.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    DataVector::DataVector(size_t size) :
      size(size), unused(0), inc_elems(100) {
      // create new vector
      this->data = new float_t[size];
    }

    DataVector::DataVector(size_t size, float_t value) :
      DataVector(size) {
      setAll(value);
    }

    DataVector::DataVector(const DataVector& vec) :
      DataVector(vec.size) {
      // copy data
      std::memcpy(this->data, vec.data, size * sizeof(float_t));
    }

    DataVector::DataVector(float_t* input, size_t size) :
      DataVector(size) {
      // copy data
      std::memcpy(this->data, input, size * sizeof(float_t));
    }

    DataVector::DataVector(std::vector<float_t> input) :
      DataVector(input.size()) {
      // copy data
      std::copy(input.begin(), input.end(), this->data);
    }

    DataVector::DataVector(std::vector<int> input) :
      DataVector(input.size()) {
      // copy data
      int in = 0;

      for (std::vector<int>::iterator it = input.begin(); it < input.end();
           it++) {
        data[in] = static_cast<float_t>(*it);
        in++;
      }
    }

    DataVector::DataVector(DataVectorDefinition& DataVectorDef) {
      setDataVectorDefinition(DataVectorDef);
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
        throw SGPP::base::algorithm_exception("more indices than entries!");
      }

      float_t* newdata = new float_t[remainingIndex.size()];

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
      float_t* newdata = new float_t[size];
      // copy entries of old vector
      std::memcpy(newdata, this->data, std::min(this->size, size)
                  * sizeof(float_t));
      delete[] this->data;

      this->data = newdata;
      this->size = size;
      this->unused = 0;
    }

    void DataVector::resizeZero(size_t size) {
      // don't do anyhing, if vector already has the correct size
      if (size == this->size) {
        return;
      }

      // create new vector
      float_t* newdata = new float_t[size];
      // copy entries of old vector
      std::memcpy(newdata, this->data, std::min(this->size, size)
                  * sizeof(float_t));

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
      float_t* newdata = new float_t[(size + add)];
      // copy entries of old vector
      std::memcpy(newdata, this->data, this->size * sizeof(float_t));

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

    size_t DataVector::append(float_t value) {
      size_t x = append();
      data[x] = value;
      return x;
    }

    void DataVector::insert(size_t index, float_t value) {
      if (index > size) {
        throw new SGPP::base::data_exception(
          "DataVector::insert : index out of bounds");
      }

      append();
      std::memmove(data + index + 1, data + index,
                   (size - 1 - index) * sizeof(float_t));
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

    void DataVector::setAll(float_t value) {
      for (size_t i = 0; i < size; i++) {
        data[i] = value;
      }
    }

    void DataVector::set(size_t i, float_t value) {
      data[i] = value;
    }

    void DataVector::copyFrom(const DataVector& vec) {
      // don't copy from yourself
      if (this == &vec) {
        return;
      }

      /*
       if (size != vec.size) {
       delete[] data;
       size = vec.size;
       this->data = new float_t[size];
       }
       */
      std::memcpy(this->data, vec.data, std::min(size, vec.size)
                  * sizeof(float_t));
    }

    DataVector& DataVector::operator=(const DataVector& vec) {
      if (this == &vec) {
        return *this;
      }

      if (size != vec.size) {
        throw new SGPP::base::data_exception(
          "DataVector::add : Dimensions do not match");
        //        delete[] data;
        //        size = vec.size;
        //        this->data = new float_t[size];
      }

      std::memcpy(this->data, vec.data, size * sizeof(float_t));
      return *this;
    }

    void DataVector::add(DataVector& vec) {
      if (size != vec.size) {
        throw new SGPP::base::data_exception(
          "DataVector::add : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] += vec.data[i];
      }
    }

    void DataVector::sub(const DataVector& vec) {
      if (size != vec.size) {
        throw new SGPP::base::data_exception(
          "DataVector::sub : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] -= vec.data[i];
      }
    }

    void DataVector::componentwise_mult(DataVector& vec) {
      if (size != vec.size) {
        throw new SGPP::base::data_exception(
          "DataVector::componentwise_mult : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] *= vec.data[i];
      }
    }

    void DataVector::componentwise_div(DataVector& vec) {
      if (size != vec.size) {
        throw new SGPP::base::data_exception(
          "DataVector::componentwise_div : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] /= vec.data[i];
      }
    }

    float_t DataVector::dotProduct(const DataVector& vec) const {
      float_t sum = 0.0;

      for (size_t i = 0; i < size; i++) {
        sum += data[i] * vec.data[i];
      }

      return sum;
    }

    void DataVector::mult(float_t scalar) {
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

    float_t DataVector::sum() const {
      float_t result = 0.0;

      for (size_t i = 0; i < size; i++) {
        result += data[i];
      }

      return result;
    }

    float_t DataVector::maxNorm() const {
      float_t max = 0.0;

      for (size_t i = 0; i < size; i++) {
        if (max < fabs(data[i])) {
          max = fabs(data[i]);
        }
      }

      return max;
    }

    float_t DataVector::RMSNorm() const {
      float_t rmsNorm;
      DataVector temp(*this);

      temp.sqr();
      rmsNorm = temp.sum();
      rmsNorm /= static_cast<float_t>(temp.getSize());
      rmsNorm = std::sqrt(rmsNorm);

      return rmsNorm;
    }

    float_t DataVector::l2Norm() const {
      float_t l2Norm;
      DataVector temp(*this);

      temp.componentwise_mult(temp);
      l2Norm = temp.sum();
      l2Norm = std::sqrt(l2Norm);

      return l2Norm;
    }

    void DataVector::partitionClasses(float_t threshold) {
      for (size_t i = 0; i < size; i++) {
        data[i] = data[i] > threshold ? 1.0 : -1.0;
      }
    }

    void DataVector::axpy(float_t a, DataVector& x) {
      if (size != x.size) {
        return;
      }

      float_t* p_x = x.data;
      float_t* p_d = data;

      for (size_t i = 0; i < size; i++) {
        p_d[i] += a * p_x[i];
      }
    }

    void DataVector::normalize() {
      normalize(0.0);
    }

    void DataVector::normalize(float_t border) {
      float_t min, max;
      minmax(&min, &max);

      float_t delta = (max - min) / (1 - 2 * border);

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

    float_t DataVector::min() const {
      float_t min = INFINITY;

      for (size_t i = 0; i < size; i++) {
        if (min > data[i]) {
          min = data[i];
        }
      }

      return min;
    }

    float_t DataVector::max() const {
      float_t max = -INFINITY;

      for (size_t i = 0; i < size; i++) {
        if (max < data[i]) {
          max = data[i];
        }
      }

      return max;
    }

    void DataVector::minmax(float_t* min, float_t* max) const {
      float_t min_t = INFINITY;
      float_t max_t = -INFINITY;

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

    float_t* DataVector::getPointer() {
      return data;
    }

    const float_t* DataVector::getPointer() const {
      return data;
    }

    DataVector::~DataVector() {
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
  }
}
