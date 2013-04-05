/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author  Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/datatypes/DataVectorSP.hpp"
#include "base/exception/data_exception.hpp"
#include "base/exception/algorithm_exception.hpp"
#include <string.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include "base/tools/AlignedMemory.hpp"
namespace sg {
  namespace base {

    DataVectorSP::DataVectorSP(size_t size) :
      size(size), unused(0), inc_elems(100) {
      // create new vector
      this->data = new float[size];
    }

    DataVectorSP::DataVectorSP(const DataVectorSP& vec) :
      unused(0), inc_elems(100) {
      this->size = vec.size;
      // create new vector
      this->data = new float[size];
      // copy data
      memcpy(this->data, vec.data, size * sizeof(float));
    }

    DataVectorSP::DataVectorSP(float* input, size_t size) :
      unused(0), inc_elems(100) {
      this->size = size;
      // create new vector
      this->data = new float[size];
      // copy data
      memcpy(this->data, input, size * sizeof(float));
    }

    void DataVectorSP::restructure(std::vector<size_t>& remainingIndex) {
      if (remainingIndex.size() > this->size) {
        throw sg::base::algorithm_exception("more indices than entries!");
      }

      float* newdata = new float[remainingIndex.size()];

      for (size_t i = 0; i < remainingIndex.size(); i++) {
        newdata[i] = this->data[remainingIndex[i]];
      }

      delete[] this->data;

      this->data = newdata;
      this->size = remainingIndex.size();
    }

    void DataVectorSP::resize(size_t size) {
      // don't do anyhing, if vector already has the correct size
      if (size == this->size) {
        return;
      }

      // create new vector
      float* newdata = new float[size];
      // copy entries of old vector
      memcpy(newdata, this->data, std::min(this->size, size) * sizeof(float));

      delete[] this->data;

      this->data = newdata;
      this->size = size;
    }

    void DataVectorSP::resizeZero(size_t size) {
      // don't do anyhing, if vector already has the correct size
      if (size == this->size) {
        return;
      }

      // create new vector
      float* newdata = new float[size];
      // copy entries of old vector
      memcpy(newdata, this->data, std::min(this->size, size) * sizeof(float));

      // set new elements to zero
      for (size_t i = std::min(this->size, size); i < size; i++) {
        newdata[i] = 0.0f;
      }

      delete[] this->data;

      this->data = newdata;
      this->size = size;
    }

    void DataVectorSP::addSize(size_t add) {
      // create new vector
      float* newdata = new float[(size + add)];
      // copy entries of old vector
      memcpy(newdata, this->data, this->size * sizeof(float));

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

      memcpy(this->data, vec.data, std::min(size, vec.size) * sizeof(float));
    }

    DataVectorSP& DataVectorSP::operator=(const DataVectorSP& vec) {
      if (this == &vec) {
        return *this;
      }

      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVectorSP::add : Dimensions do not match");
      }

      memcpy(this->data, vec.data, size * sizeof(float));
      return *this;
    }

    void DataVectorSP::add(DataVectorSP& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVectorSP::add : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] += vec.data[i];
      }
    }

    void DataVectorSP::sub(const DataVectorSP& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVectorSP::sub : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] -= vec.data[i];
      }
    }

    void DataVectorSP::componentwise_mult(DataVectorSP& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVectorSP::componentwise_mult : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] *= vec.data[i];
      }
    }

    void DataVectorSP::componentwise_div(DataVectorSP& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVectorSP::componentwise_div : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] /= vec.data[i];
      }
    }

    float DataVectorSP::dotProduct(DataVectorSP& vec) {
      float sum = 0.0;

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

    float DataVectorSP::sum() {
      float result = 0.0f;

      for (size_t i = 0; i < size; i++) {
        result += data[i];
      }

      return result;
    }

    float DataVectorSP::maxNorm() {
      float max = 0.0f;

      for (size_t i = 0; i < size; i++) {
        if (max < static_cast<float>(fabs(data[i]))) {
          max = static_cast<float>(fabs(data[i]));
        }
      }

      return max;
    }

    float DataVectorSP::RMSNorm() {
      float l2Norm;
      DataVectorSP temp(*this);

      temp.componentwise_mult(temp);
      l2Norm = temp.sum();
      l2Norm /= static_cast<float>(temp.getSize());
      l2Norm = std::sqrt(l2Norm);

      return l2Norm;
    }

    float DataVectorSP::l2Norm() {
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

    void DataVectorSP::toString(std::string& text) {
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

    std::string DataVectorSP::toString() {
      std::string str;
      toString(str);
      return str;
    }

    float DataVectorSP::min() {
      float min = data[0];

      for (size_t i = 1; i < size; i++) {
        if (min > data[i]) {
          min = data[i];
        }
      }

      return min;
    }

    float DataVectorSP::max() {
      float max = data[0];

      for (size_t i = 1; i < size; i++) {
        if (max < data[i]) {
          max = data[i];
        }
      }

      return max;
    }

    void DataVectorSP::minmax(float* min, float* max) {
      float min_t = data[0];
      float max_t = data[0];

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

    DataVectorSP::~DataVectorSP() {
      delete[] data;
    }

    size_t DataVectorSP::getNumberNonZero() {
      size_t nonZero = 0;

      for (size_t i = 0; i < size; i++) {
        if (fabs(data[i]) > 0.0f) {
          nonZero++;
        }
      }

      return nonZero;
    }
  }
}
