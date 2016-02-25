// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

namespace SGPP {
namespace base {

DataMatrix::DataMatrix() : nrows(0), ncols(0), unused(0), inc_rows(100) {
  // create new vector
  this->data = new float_t[0];
}

DataMatrix::DataMatrix(size_t nrows, size_t ncols)
    : nrows(nrows), ncols(ncols), unused(0), inc_rows(100) {
  // create new vector
  this->data = new float_t[nrows * ncols];
}

DataMatrix::DataMatrix(size_t nrows, size_t ncols, float_t value) : DataMatrix(nrows, ncols) {
  setAll(value);
}

DataMatrix::DataMatrix(const DataMatrix& matr) : DataMatrix(matr.nrows, matr.ncols) {
  // copy data
  std::memcpy(this->data, matr.data, nrows * ncols * sizeof(float_t));
}

DataMatrix::DataMatrix(float_t* input, size_t nrows, size_t ncols) : DataMatrix(nrows, ncols) {
  // copy data
  std::memcpy(this->data, input, nrows * ncols * sizeof(float_t));
}

DataMatrix DataMatrix::fromFile(const std::string& fileName) {
  std::ifstream f(fileName, std::ifstream::in);
  f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  std::string content;
  content.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
  return DataMatrix::fromString(content);
}

DataMatrix DataMatrix::fromString(const std::string& serializedVector) {
  DataMatrix m;

  enum class PARSER_STATE { INIT, ROW, ROWVALUE, ROWCOMMAEND, COMMAEND, END };

  PARSER_STATE state = PARSER_STATE::INIT;

  DataVector row;

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
      state = PARSER_STATE::ROW;
      i += 1;
    } else if (state == PARSER_STATE::ROW) {
      row.resize(0);
      state = PARSER_STATE::ROWVALUE;
      i++;
    } else if (state == PARSER_STATE::ROWVALUE) {
//      size_t next;
#if USE_DOUBLE_PRECISION == 1
      //      double value = std::stod(&(serializedVector[i]), &next);
      double value = std::atof(&(serializedVector[i]));
#else
      float value = std::stof(&(serializedVector[i]));
#endif
      row.append(value);
      state = PARSER_STATE::ROWCOMMAEND;
      //      i += next;
      while (serializedVector[i] != ',' && serializedVector[i] != ']') i++;
    } else if (state == PARSER_STATE::ROWCOMMAEND) {
      if (c == ',') {
        state = PARSER_STATE::ROWVALUE;
        i++;
      } else if (c == ']') {
        if (m.getNrows() == 0) {
          // set up the dimension after having read the first row
          m.resize(0, row.getSize());
        }
        m.appendRow(row);
        state = PARSER_STATE::COMMAEND;
        i++;
      }
    } else if (state == PARSER_STATE::COMMAEND) {
      if (c == ',') {
        state = PARSER_STATE::ROW;
        i++;
      } else if (c == ']') {
        state = PARSER_STATE::END;
        i++;
      }
    } else if (state == PARSER_STATE::END) {
      // only reached if a non-whitespace character was encountered after closing brace
      throw data_exception("error: could not parse DataMatrix file");
    }
  }
  return m;
}

/**
 DataMatrix::DataMatrix(DataMatrixDefinition& DataMatrixDef) {
 this->nrows = DataMatrixDef.nrows;
 this->ncols = DataMatrixDef.ncols;
 this->unused = DataMatrixDef.unused;
 this->data = DataMatrixDef.pointerToData;
 }

 void DataMatrix::getDataMatrixDefinition(DataMatrixDefinition& DataMatrixDef) {
 DataMatrixDef.nrows = this->nrows;
 DataMatrixDef.ncols = this->ncols;
 DataMatrixDef.unused = this->unused;
 DataMatrixDef.pointerToData = this->data;
 }
 **/

void DataMatrix::resize(size_t nrows) {
  // don't do anything, if matrix already has the correct size
  if (nrows == this->nrows) {
    return;
  }

  // create new matrix
  float_t* newdata = new float_t[nrows * this->ncols];
  // copy entries of old matrix
  std::memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols * sizeof(float_t));
  delete[] this->data;

  this->data = newdata;
  this->nrows = nrows;
  this->unused = 0;
}

void DataMatrix::resize(size_t nrows, size_t ncols) {
  // don't do anything, if matrix already has the correct size
  if ((nrows == this->nrows) && (ncols == this->ncols)) {
    return;
  }

  // don't copy data, if matrix already has the correct number of entries
  if (this->nrows * this->ncols != nrows * ncols) {
    // create new matrix
    float_t* newdata = new float_t[nrows * ncols];
    // copy entries of old matrix
    std::memcpy(newdata, this->data,
                std::min(this->nrows * this->ncols, nrows * ncols) * sizeof(float_t));
    delete[] this->data;
    this->data = newdata;
  }

  this->nrows = nrows;
  this->ncols = ncols;
  this->unused = 0;
}

void DataMatrix::resizeZero(size_t nrows) {
  // don't do anything, if matrix already has the correct size
  if (nrows == this->nrows) {
    return;
  }

  // create new matrix
  float_t* newdata = new float_t[nrows * this->ncols];
  // copy entries of old matrix
  std::memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols * sizeof(float_t));

  // set new elements to zero
  for (size_t i = std::min(this->nrows, nrows) * this->ncols; i < nrows * this->ncols; i++) {
    newdata[i] = 0.0;
  }

  delete[] this->data;

  this->data = newdata;
  this->nrows = nrows;
  this->unused = 0;
}

void DataMatrix::resizeZero(size_t nrows, size_t ncols) {
  // don't do anything, if matrix already has the correct size
  if ((nrows == this->nrows) && (ncols == this->ncols)) {
    return;
  }

  // don't copy data, if matrix already has the correct number of entries
  if (this->nrows * this->ncols != nrows * ncols) {
    // create new matrix
    float_t* newdata = new float_t[nrows * ncols];
    // copy entries of old matrix
    std::memcpy(newdata, this->data,
                std::min(this->nrows * this->ncols, nrows * ncols) * sizeof(float_t));

    // set new elements to zero
    for (size_t i = std::min(this->nrows * this->ncols, nrows * ncols); i < nrows * ncols; i++) {
      newdata[i] = 0.0;
    }

    delete[] this->data;
    this->data = newdata;
  }

  this->nrows = nrows;
  this->ncols = ncols;
  this->unused = 0;
}

void DataMatrix::addSize(size_t inc_nrows) {
  // create new vector
  float_t* newdata = new float_t[(this->nrows + inc_nrows) * this->ncols];
  // copy entries of old vector
  std::memcpy(newdata, this->data, this->nrows * this->ncols * sizeof(float_t));

  delete[] this->data;

  this->data = newdata;
  this->unused = inc_nrows;
}

size_t DataMatrix::appendRow() {
  // enlarge, if necessary
  if (unused == 0) {
    addSize(this->inc_rows);
  }

  size_t x = nrows;
  nrows++;
  unused--;

  return x;
}

void DataMatrix::transpose() {
  float_t* newData = new float_t[nrows * ncols];

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
  unused = 0;
}

size_t DataMatrix::appendRow(const DataVector& vec) {
  if (vec.getSize() != this->ncols) {
    throw SGPP::base::data_exception("DataMatrix::appendRow : Dimensions do not match");
  }

  size_t x = appendRow();
  // copy data
  std::memcpy(&this->data[x * this->ncols], vec.getPointer(), this->ncols * sizeof(float_t));
  return x;
}

void DataMatrix::setAll(float_t value) {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] = value;
  }
}

void DataMatrix::getRow(size_t row, DataVector& vec) const {
  if (vec.getSize() != this->ncols) {
    throw SGPP::base::data_exception("DataMatrix::getRow : Dimensions do not match");
  }

  for (size_t i = 0; i < this->ncols; i++) {
    vec[i] = this->data[row * ncols + i];
  }
}

void DataMatrix::getRow(size_t row, std::vector<float_t>& vec) const {
  vec.clear();

  for (size_t i = 0; i < this->ncols; i++) {
    vec.push_back(data[row * ncols + i]);
  }
}

void DataMatrix::setRow(size_t row, const DataVector& vec) {
  if (vec.getSize() != this->ncols) {
    throw SGPP::base::data_exception("DataMatrix::setRow : Dimensions do not match");
  }

  for (size_t i = 0; i < this->ncols; i++) {
    this->data[row * ncols + i] = vec.get(i);
  }
}

void DataMatrix::getColumn(size_t col, DataVector& vec) const {
  if (vec.getSize() != this->nrows) {
    throw SGPP::base::data_exception("DataMatrix::getColumn : Dimensions do not match");
  }

  for (size_t j = 0; j < this->nrows; j++) {
    vec[j] = data[j * ncols + col];
  }
}

void DataMatrix::setColumn(size_t col, const DataVector& vec) {
  if (vec.getSize() != this->nrows) {
    throw SGPP::base::data_exception("DataMatrix::setColumn : Dimensions do not match");
  }

  for (size_t j = 0; j < this->nrows; j++) {
    data[j * ncols + col] = vec.get(j);
  }
}

void DataMatrix::copyFrom(const DataMatrix& matr) {
  // don't copy from yourself
  if (this == &matr) {
    return;
  }

  /*
   if (nrows != vec.nrows || ncols != vec.ncols) {
   delete[] data;
   nrows = vec.nrows;
   ncols = vec.ncols;
   this->data = new float_t[nrows * ncols];
   }
   */
  std::memcpy(this->data, matr.data,
              std::min(this->nrows * this->ncols, matr.nrows * matr.ncols) * sizeof(float_t));
}

/*
 void DataMatrix::copySmall(const DataMatrix& vec) {
 if (this == &vec) {
 return;
 }

 if (vec.ncols != 1 || ncols != 1 || nrows < vec.nrows) {
 return;
 }
 std::memcpy(this->data, vec.data, vec.nrows * sizeof(float_t));
 }

 DataMatrix& DataMatrix::operator=(const DataMatrix &vec) {
 if (this == &vec) {
 return *this;
 }

 if (nrows != vec.nrows || ncols != vec.ncols) {
 delete[] data;
 nrows = vec.nrows;
 ncols = vec.ncols;
 this->data = new float_t[nrows * ncols];
 }
 std::memcpy(this->data, vec.data, nrows * ncols * sizeof(float_t));
 return *this;
 }
 */
DataMatrix& DataMatrix::operator=(const DataMatrix& matr) {
  if (this == &matr) {
    return *this;
  }

  if (nrows * ncols != matr.ncols * matr.nrows) {
    throw SGPP::base::data_exception("DataMatrix::= : Dimensions do not match");
  }

  std::memcpy(this->data, matr.data, nrows * ncols * sizeof(float_t));
  return *this;
}

void DataMatrix::add(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw SGPP::base::data_exception("DataMatrix::add : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] += matr.data[i];
  }
}

void DataMatrix::sub(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw SGPP::base::data_exception("DataMatrix::sub : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] -= matr.data[i];
  }
}

void DataMatrix::addReduce(DataVector& reduction) {
  if (this->nrows != reduction.getSize()) {
    throw SGPP::base::data_exception("DataMatrix::addReduce : Dimensions do not match");
  }

  for (size_t i = 0; i < this->nrows; i++) {
    float_t tmp = 0.0;

    for (size_t j = 0; j < this->ncols; j++) {
      tmp += this->data[(i * this->ncols) + j];
    }

    reduction.set(i, tmp);
  }
}

void DataMatrix::addReduce(DataVector& reduction, DataVector& beta, size_t start_beta) {
  if (this->nrows != reduction.getSize()) {
    throw SGPP::base::data_exception("DataMatrix::addReduce : Dimensions do not match (reduction)");
  }

  if (this->ncols + start_beta > beta.getSize()) {
    throw SGPP::base::data_exception("DataMatrix::addReduce : Dimensions do not match (beta)");
  }

  for (size_t i = 0; i < this->nrows; i++) {
    float_t tmp = 0.0;

    for (size_t j = 0; j < this->ncols; j++) {
      tmp += beta[j + start_beta] * this->data[(i * this->ncols) + j];
    }

    reduction.set(i, reduction[i] + tmp);
  }
}

void DataMatrix::expand(const DataVector& expand) {
  if (this->nrows != expand.getSize()) {
    throw SGPP::base::data_exception("DataMatrix::expand : Dimensions do not match");
  }

  for (size_t i = 0; i < this->nrows; i++) {
    for (size_t j = 0; j < this->ncols; j++) {
      this->data[(i * this->ncols) + j] = expand.get(i);
    }
  }
}

void DataMatrix::componentwise_mult(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw SGPP::base::data_exception("DataMatrix::componentwise_mult : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] *= matr.data[i];
  }
}

void DataMatrix::componentwise_div(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw SGPP::base::data_exception("DataMatrix::componentwise_div : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] /= matr.data[i];
  }
}

/*
 void DataMatrix::getLine(int row, DataMatrix& vec) {
 for (int i = 0; i < this->ncols; i++) {
 vec[i] = data[row * ncols + i];
 }
 }

 void DataMatrix::getLine(int row, std::vector<float_t>& vec) {
 vec.clear();

 for (int i = 0; i < this->ncols; i++) {
 vec.push_back(data[row * ncols + i]);
 }
 }
 */

/*
 float_t DataMatrix::dotProduct(DataMatrix &vec) {
 float_t sum = 0.0;

 for (int i = 0; i < nrows; i++) {
 sum += data[i] * vec.data[i];
 }
 return sum;
 }
 */

void DataMatrix::mult(float_t scalar) {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] *= scalar;
  }
}

void DataMatrix::mult(const DataVector& x, DataVector& y) {
  if (ncols != x.getSize()) {
    throw SGPP::base::data_exception("DataMatrix::mult : Dimensions do not match (x)");
  }

  if (nrows != y.getSize()) {
    throw SGPP::base::data_exception("DataMatrix::mult : Dimensions do not match (y)");
  }

  for (size_t i = 0; i < nrows; i++) {
    float_t entry = 0.0;

    for (size_t j = 0; j < ncols; j++) {
      entry += data[(i * ncols) + j] * x[j];
    }

    y[i] = entry;
  }
}

void DataMatrix::sqr() {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] = data[i] * data[i];
  }
}

void DataMatrix::sqrt() {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] = std::sqrt(data[i]);
  }
}

void DataMatrix::abs() {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] = std::fabs(data[i]);
  }
}

float_t DataMatrix::sum() const {
  size_t n = nrows * ncols;
  float_t result = 0.0;

  for (size_t i = 0; i < n; i++) {
    result += data[i];
  }

  return result;
}
/*
 float_t DataMatrix::maxNorm() {
 int n = nrows * ncols;
 float_t max = 0.0;
 for (int i = 0; i < n; i++) {
 if (max < fabs(data[i]))
 {
 max = fabs(data[i]);
 }
 }
 return max;
 }

 void DataMatrix::partitionClasses(float_t border) {
 int n = nrows * ncols;
 for (int i = 0; i < n; i++) {
 data[i] = data[i] > border ? 1.0 : -1.0;
 }
 }

 void DataMatrix::axpy(float_t alpha, DataMatrix& x) {
 if (nrows != x.nrows || ncols != x.ncols) {
 return;
 }
 int n = nrows * ncols;
 float_t* p_x = x.data;
 float_t* p_d = data;

 for (int i = 0; i < n; i++) {
 p_d[i] += alpha * p_x[i];
 }
 }
 */

void DataMatrix::normalizeDimension(size_t d) { normalizeDimension(d, 0.0); }

void DataMatrix::normalizeDimension(size_t d, float_t border) {
  size_t n = nrows * ncols;

  if (ncols <= d) {
    throw SGPP::base::data_exception(
        "DataMatrix::normalizeDimension : Not enough columns in DataMatrix");
  }

  // determine min and max
  float_t xmin, xmax;
  minmax(d, &xmin, &xmax);

  float_t delta = (xmax - xmin) / (1 - 2 * border);

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

void DataMatrix::toString(std::string& text) const {
  std::stringstream str;
  str << "[";

  for (size_t i = 0; i < nrows; i++) {
    str << "[";

    for (size_t j = 0; j < ncols; j++) {
      if (j != 0) {
        str << ", ";
      }

      str << data[i * ncols + j];
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

std::string DataMatrix::toString() const {
  std::string str;
  toString(str);
  return str;
}

void DataMatrix::toFile(const std::string& fileName) const {
  std::ofstream f(fileName, std::ofstream::out);
  f << this->toString();
  f.close();
}

float_t DataMatrix::min(size_t d) const {
  size_t n = nrows * ncols;
  float_t min = INFINITY;

  for (size_t i = d; i < n; i += ncols) {
    if (min > data[i]) {
      min = data[i];
    }
  }

  return min;
}

float_t DataMatrix::min() const {
  size_t n = nrows * ncols;
  float_t min = INFINITY;

  for (size_t i = 0; i < n; i++) {
    if (min > data[i]) {
      min = data[i];
    }
  }

  return min;
}

float_t DataMatrix::max(size_t d) const {
  size_t n = nrows * ncols;
  float_t max = -INFINITY;

  for (size_t i = d; i < n; i += ncols) {
    if (max < data[i]) {
      max = data[i];
    }
  }

  return max;
}

float_t DataMatrix::max() const {
  size_t n = nrows * ncols;
  float_t max = -INFINITY;

  for (size_t i = 0; i < n; i++) {
    if (max < data[i]) {
      max = data[i];
    }
  }

  return max;
}

void DataMatrix::minmax(size_t col, float_t* min, float_t* max) const {
  size_t n = nrows * ncols;

  if (ncols <= col) {
    throw SGPP::base::data_exception("DataMatrix::minmax : Not enough entries in DataMatrix");
  }

  // find min and max of column col
  float_t min_t = INFINITY;
  float_t max_t = -INFINITY;

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

void DataMatrix::minmax(float_t* min, float_t* max) const {
  size_t n = nrows * ncols;

  float_t min_t = INFINITY;
  float_t max_t = -INFINITY;

  for (size_t i = 0; i < n; i++) {
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

float_t* DataMatrix::getPointer() { return data; }

const float_t* DataMatrix::getPointer() const { return data; }

DataMatrix::~DataMatrix() { delete[] data; }

size_t DataMatrix::getNumberNonZero() const {
  size_t n = nrows * ncols;
  size_t nonZero = 0;

  for (size_t i = 0; i < n; i++) {
    if (fabs(data[i]) > 0.0) {
      nonZero++;
    }
  }

  return nonZero;
}

}  // namespace base
}  // namespace SGPP
