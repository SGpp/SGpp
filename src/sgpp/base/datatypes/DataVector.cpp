/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de)

#include "base/datatypes/DataVector.hpp"
#include "base/exception/data_exception.hpp"
#include "base/exception/algorithm_exception.hpp"
#include <string.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include "base/tools/AlignedMemory.hpp"

#include <iostream>

namespace sg {
  namespace base {

    DataVector::DataVector(size_t size) :
      size(size), unused(0), inc_elems(100) {
      // create new vector
      this->data = new double[size];
    }

    DataVector::DataVector(const DataVector& vec) :
      unused(0), inc_elems(100) {
      this->size = vec.size;
      // create new vector
      this->data = new double[size];
      // copy data
      memcpy(this->data, vec.data, size * sizeof(double));
    }

    DataVector::DataVector(double* input, size_t size) :
      unused(0), inc_elems(100) {
      this->size = size;
      // create new vector
      this->data = new double[size];
      // copy data
      memcpy(this->data, input, size * sizeof(double));
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
        throw sg::base::algorithm_exception("more indices than entries!");
      }

      double* newdata = new double[remainingIndex.size()];

      for (size_t i = 0; i < remainingIndex.size(); i++) {
        newdata[i] = this->data[remainingIndex[i]];
      }

      delete[] this->data;

      this->data = newdata;
      this->size = remainingIndex.size();
    }

    void DataVector::resize(size_t size) {
      // don't do anyhing, if vector already has the correct size
      if (size == this->size) {
        return;
      }

      // create new vector
      double* newdata = new double[size];
      // copy entries of old vector
      memcpy(newdata, this->data, std::min(this->size, size) * sizeof(double));
      delete[] this->data;

      this->data = newdata;
      this->size = size;
    }

    void DataVector::resizeZero(size_t size) {
      // don't do anyhing, if vector already has the correct size
      if (size == this->size) {
        return;
      }

      // create new vector
      double* newdata = new double[size];
      // copy entries of old vector
      memcpy(newdata, this->data, std::min(this->size, size) * sizeof(double));

      // set new elements to zero
      for (size_t i = std::min(this->size, size); i < size; i++) {
        newdata[i] = 0.0;
      }

      delete[] this->data;

      this->data = newdata;
      this->size = size;
    }

    void DataVector::addSize(size_t add) {
      // create new vector
      double* newdata = new double[(size + add)];
      // copy entries of old vector
      memcpy(newdata, this->data, this->size * sizeof(double));

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

    void DataVector::set(size_t i, double value) {
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
          this->data = new double[size];
      }
      */
      memcpy(this->data, vec.data, std::min(size, vec.size) * sizeof(double));
    }

    DataVector& DataVector::operator=(const DataVector& vec) {
      if (this == &vec) {
        return *this;
      }

      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVector::add : Dimensions do not match");
        //        delete[] data;
        //        size = vec.size;
        //        this->data = new double[size];
      }

      memcpy(this->data, vec.data, size * sizeof(double));
      return *this;
    }

    void DataVector::add(DataVector& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVector::add : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] += vec.data[i];
      }
    }

    void DataVector::sub(const DataVector& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVector::sub : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] -= vec.data[i];
      }
    }

    void DataVector::componentwise_mult(DataVector& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVector::componentwise_mult : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] *= vec.data[i];
      }
    }

    void DataVector::componentwise_div(DataVector& vec) {
      if (size != vec.size) {
        throw new sg::base::data_exception(
          "DataVector::componentwise_div : Dimensions do not match");
      }

      for (size_t i = 0; i < size; i++) {
        data[i] /= vec.data[i];
      }
    }

    double DataVector::dotProduct(DataVector& vec) {
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

    double DataVector::sum() {
      double result = 0.0;

      for (size_t i = 0; i < size; i++) {
        result += data[i];
      }

      return result;
    }

    double DataVector::maxNorm() {
      double max = 0.0;

      for (size_t i = 0; i < size; i++) {
        if (max < fabs(data[i])) {
          max = fabs(data[i]);
        }
      }

      return max;
    }

    double DataVector::RMSNorm() {
      double rmsNorm;
      DataVector temp(*this);

      temp.sqr();
      rmsNorm = temp.sum();
      rmsNorm /= static_cast<double>(temp.getSize());
      rmsNorm = std::sqrt(rmsNorm);

      return rmsNorm;
    }

    double DataVector::l2Norm() {
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

    void DataVector::normalize() {
      normalize(0.0);
    }

    void DataVector::normalize(double border) {
      double min, max;
      minmax(&min, &max);

      double delta = (max - min) / (1 - 2 * border);

      for (size_t i = 0; i < size; i++) {
        data[i] = (data[i] - min) / delta + border;
      }
    }

    void DataVector::toString(std::string& text) {
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

    std::string DataVector::toString() {
      std::string str;
      toString(str);
      return str;
    }

    double DataVector::min() {
      double min = data[0];

      for (size_t i = 1; i < size; i++) {
        if (min > data[i]) {
          min = data[i];
        }
      }

      return min;
    }

    double DataVector::max() {
      double max = data[0];

      for (size_t i = 1; i < size; i++) {
        if (max < data[i]) {
          max = data[i];
        }
      }

      return max;
    }

    void DataVector::minmax(double* min, double* max) {
      double min_t = data[0];
      double max_t = data[0];

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

    double* DataVector::getPointer() {
      return data;
    }

    DataVector::~DataVector() {
      delete[] data;
    }

    size_t DataVector::getNumberNonZero() {
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
