/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/datatypes/DataVectorSP.hpp"
#include "base/datatypes/DataMatrixSP.hpp"
#include "base/exception/data_exception.hpp"

#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>

#include <iostream>

#include "base/tools/AlignedMemory.hpp"
namespace sg {
  namespace base {

    DataMatrixSP::DataMatrixSP(size_t nrows, size_t ncols) :
      nrows(nrows), ncols(ncols), unused(0), inc_rows(100) {
      // create new vector
      this->data = new float[nrows * ncols];
    }

    DataMatrixSP::DataMatrixSP(const DataMatrixSP& matr) :
      unused(0), inc_rows(100) {
      this->nrows = matr.nrows;
      this->ncols = matr.ncols;
      // create new vector
      this->data = new float[nrows * ncols];
      // copy data
      memcpy(this->data, matr.data, nrows * ncols * sizeof(float));
    }

    DataMatrixSP::DataMatrixSP(float* input, size_t nrows, size_t ncols) :
      nrows(nrows), ncols(ncols), unused(0), inc_rows(100) {
      // create new vector
      this->data = new float[nrows * ncols];
      // copy data
      memcpy(this->data, input, nrows * ncols * sizeof(float));
    }

    void DataMatrixSP::resize(size_t nrows) {
      // don't do anyhing, if vector already has the correct size
      if (nrows == this->nrows) {
        return;
      }

      // create new vector
      float* newdata = new float[nrows * this->ncols];
      // copy entries of old vector
      memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols
             * sizeof(float));
      delete[] this->data;

      this->data = newdata;
      this->nrows = nrows;
    }

    void DataMatrixSP::resizeZero(size_t nrows) {
      // don't do anyhing, if vector already has the correct size
      if (nrows == this->nrows) {
        return;
      }

      // create new vector
      float* newdata = new float[nrows * this->ncols];
      // copy entries of old vector
      memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols
             * sizeof(float));

      // set new elements to zero
      for (size_t i = std::min(this->nrows, nrows) * this->ncols; i < nrows
           * this->ncols; i++) {
        newdata[i] = 0.0;
      }

      delete[] this->data;

      this->data = newdata;
      this->nrows = nrows;
    }

    void DataMatrixSP::addSize(size_t inc_rows) {
      // create new vector
      float* newdata = new float[(this->nrows + inc_rows) * this->ncols];
      // copy entries of old vector
      memcpy(newdata, this->data, this->nrows * this->ncols * sizeof(float));

      delete[] this->data;

      this->data = newdata;
      this->unused = inc_rows;
    }

    size_t DataMatrixSP::appendRow() {
      // enlarge, if necessary
      if (unused == 0) {
        addSize(this->inc_rows);
      }

      size_t x = nrows;
      nrows++;
      unused--;

      return x;
    }

    void DataMatrixSP::transpose() {
      float* newData = new float[nrows * ncols];

      for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
          newData[(j * nrows) + i] = data[(i * ncols) + j];
        }
      }

      delete[] data;
      data = newData;
      size_t tmpRows = nrows;
      nrows = ncols;
      ncols = tmpRows;
    }

    size_t DataMatrixSP::appendRow(DataVectorSP& vec) {
      if (vec.getSize() != this->ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::appendRow : Dimensions do not match");
      }

      size_t x = appendRow();
      // copy data
      memcpy(&this->data[x * this->ncols], vec.getPointer(), this->ncols
             * sizeof(float));
      return x;
    }

    void DataMatrixSP::setAll(float value) {
      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] = value;
      }
    }

    void DataMatrixSP::getRow(size_t row, DataVectorSP& vec) {
      if (vec.getSize() != this->ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::getRow : Dimensions do not match");
      }

      for (size_t i = 0; i < this->ncols; i++) {
        vec[i] = this->data[row * ncols + i];
      }
    }

    void DataMatrixSP::getRow(size_t row, std::vector<float>& vec) {
      vec.clear();

      for (size_t i = 0; i < this->ncols; i++) {
        vec.push_back(data[row * ncols + i]);
      }
    }

    void DataMatrixSP::setRow(size_t row, DataVectorSP& vec) {
      if (vec.getSize() != this->ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::setRow : Dimensions do not match");
      }

      for (size_t i = 0; i < this->ncols; i++) {
        this->data[row * ncols + i] = vec[i];
      }
    }

    void DataMatrixSP::getColumn(size_t col, DataVectorSP& vec) {
      if (vec.getSize() != this->nrows) {
        throw new sg::base::data_exception(
          "DataMatrixSP::getColumn : Dimensions do not match");
      }

      for (size_t j = 0; j < this->nrows; j++) {
        vec[j] = data[j * ncols + col];
      }
    }

    void DataMatrixSP::setColumn(size_t col, DataVectorSP& vec) {
      if (vec.getSize() != this->nrows) {
        throw new sg::base::data_exception(
          "DataMatrixSP::setColumn : Dimensions do not match");
      }

      for (size_t j = 0; j < this->nrows; j++) {
        data[j * ncols + col] = vec[j];
      }
    }

    void DataMatrixSP::copyFrom(const DataMatrixSP& matr) {
      // don't copy from yourself
      if (this == &matr) {
        return;
      }

      /*
       if (nrows != vec.nrows || ncols != vec.ncols) {
       delete[] data;
       nrows = vec.nrows;
       ncols = vec.ncols;
       this->data = new float[nrows * ncols];
       }
       */
      memcpy(this->data, matr.data, std::min(this->nrows * this->ncols, matr.nrows * matr.ncols) * sizeof(float));
    }

    DataMatrixSP& DataMatrixSP::operator=(const DataMatrixSP& matr) {
      if (this == &matr) {
        return *this;
      }

      if (nrows * ncols != matr.ncols * matr.nrows) {
        throw new sg::base::data_exception(
          "DataMatrixSP::= : Dimensions do not match");
      }

      memcpy(this->data, matr.data, nrows * ncols * sizeof(float));
      return *this;
    }


    void DataMatrixSP::add(DataMatrixSP& matr) {
      if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::add : Dimensions do not match");
      }

      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] += matr.data[i];
      }
    }

    void DataMatrixSP::sub(const DataMatrixSP& matr) {
      if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::sub : Dimensions do not match");
      }

      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] -= matr.data[i];
      }
    }

    void DataMatrixSP::addReduce(DataVectorSP& reduction) {
      if (this->nrows != reduction.getSize() ) {
        throw new sg::base::data_exception(
          "DataMatrixSP::addReduce : Dimensions do not match");
      }

      for (size_t i = 0; i < this->nrows; i++) {
        float tmp = 0.0;

        for (size_t j = 0; j < this->ncols; j++) {
          tmp += this->data[(i * this->ncols) + j];
        }

        reduction.set(i, tmp);
      }
    }

    void DataMatrixSP::expand(const DataVectorSP& expand) {
      if (this->nrows != expand.getSize() ) {
        throw new sg::base::data_exception(
          "DataMatrixSP::expand : Dimensions do not match");
      }

      for (size_t i = 0; i < this->nrows; i++) {
        for (size_t j = 0; j < this->ncols; j++) {
          this->data[(i * this->ncols) + j] = expand.get(i);
        }
      }
    }

    void DataMatrixSP::componentwise_mult(DataMatrixSP& matr) {
      if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::componentwise_mult : Dimensions do not match");
      }

      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] *= matr.data[i];
      }
    }

    void DataMatrixSP::componentwise_div(DataMatrixSP& matr) {
      if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
        throw new sg::base::data_exception(
          "DataMatrixSP::componentwise_div : Dimensions do not match");
      }

      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] /= matr.data[i];
      }
    }

    void DataMatrixSP::mult(float scalar) {
      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] *= scalar;
      }
    }

    void DataMatrixSP::sqr() {
      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] = data[i] * data[i];
      }
    }

    void DataMatrixSP::sqrt() {
      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] = std::sqrt(data[i]);
      }
    }

    void DataMatrixSP::abs() {
      size_t n = nrows * ncols;

      for (size_t i = 0; i < n; i++) {
        data[i] = std::fabs(data[i]);
      }
    }

    float DataMatrixSP::sum() {
      size_t n = nrows * ncols;
      float result = 0.0;

      for (size_t i = 0; i < n; i++) {
        result += data[i];
      }

      return result;
    }

    void DataMatrixSP::normalizeDimension(size_t d) {
      normalizeDimension(d, 0.0);
    }

    void DataMatrixSP::normalizeDimension(size_t d, float border) {
      size_t n = nrows * ncols;

      if (ncols <= d) {
        throw new sg::base::data_exception(
          "DataMatrixSP::normalizeDimension : Not enough columns in DataMatrixSP");
      }

      // determine min and max
      float xmin, xmax;
      minmax(d, &xmin, &xmax);

      float delta = (xmax - xmin) / (1 - 2 * border);

      if (delta == 0.0) {
        for (size_t i = d; i < n; i += ncols) {
          data[i] = 0.5;
        }
      } else {
        for (size_t i = d; i < n; i += ncols) {
          data[i] = (data[i] - xmin) / delta + border;
        }
      }
    }

    void DataMatrixSP::toString(std::string& text) {
      std::stringstream str;
      str << "[";

      for (size_t i = 0; i < nrows; i++) {
        str << "[";

        for (size_t j = 0; j < ncols; j++) {
          if (j != 0) {
            str << ",";
          }

          str << " " << data[i * ncols + j];
        }

        if (i == nrows - 1) {
          str << "]";
        } else {
          str << "]," << std::endl;
        }
      }

      str << "]";
      text = str.str();
    }

    std::string DataMatrixSP::toString() {
      std::string str;
      toString(str);
      return str;
    }

    float DataMatrixSP::min(size_t d) {
      size_t n = nrows * ncols;
      float min = data[d];

      for (size_t i = d; i < n; i += ncols) {
        if (min > data[i]) {
          min = data[i];
        }
      }

      return min;
    }

    float DataMatrixSP::min() {
      size_t n = nrows * ncols;
      float min = data[0];

      for (size_t i = 1; i < n; i += 1) {
        if (min > data[i]) {
          min = data[i];
        }
      }

      return min;
    }

    float DataMatrixSP::max(size_t d) {
      size_t n = nrows * ncols;
      float max = data[d];

      for (size_t i = d; i < n; i += ncols) {
        if (max < data[i]) {
          max = data[i];
        }
      }

      return max;
    }

    float DataMatrixSP::max() {
      size_t n = nrows * ncols;
      float max = data[0];

      for (size_t i = 1; i < n; i += 1) {
        if (max < data[i]) {
          max = data[i];
        }
      }

      return max;
    }

    void DataMatrixSP::minmax(size_t col, float* min, float* max) {
      size_t n = nrows * ncols;

      if (ncols <= col) {
        throw new sg::base::data_exception(
          "DataMatrixSP::minmax : Not enough entries in DataMatrixSP");
      }

      // find min and max of column col
      float min_t = data[col];
      float max_t = data[col];

      for (size_t i = col; i < n; i += ncols) {
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

    void DataMatrixSP::minmax(float* min, float* max) {
      size_t n = nrows * ncols;

      if (n == 0) {
        throw new sg::base::data_exception(
          "DataMatrixSP::minmax : Empty DataMatrixSP");
      }

      float min_t = data[0];
      float max_t = data[0];

      for (size_t i = 1; i < n; i += 1) {
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

    float* DataMatrixSP::getPointer() {
      return data;
    }

    DataMatrixSP::~DataMatrixSP() {
      delete[] data;
    }

    size_t DataMatrixSP::getNumberNonZero() {
      size_t n = nrows * ncols;
      size_t nonZero = 0;

      for (size_t i = 0; i < n; i++) {
        if (fabs(data[i]) > 0.0) {
          nonZero++;
        }
      }

      return nonZero;
    }
  }
}
