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

namespace sgpp {
namespace base {

DataMatrix::DataMatrix() : nrows(0), ncols(0), unused(0), inc_rows(100) {
  // create new vector
  this->data = new double[0];
}

DataMatrix::DataMatrix(size_t nrows, size_t ncols)
    : nrows(nrows), ncols(ncols), unused(0), inc_rows(100) {
  // create new vector
  this->data = new double[nrows * ncols];
}

DataMatrix::DataMatrix(size_t nrows, size_t ncols, double value) : DataMatrix(nrows, ncols) {
  setAll(value);
}

DataMatrix::DataMatrix(const DataMatrix& matr) : DataMatrix(matr.nrows, matr.ncols) {
  // copy data
  std::memcpy(this->data, matr.data, nrows * ncols * sizeof(double));
}

DataMatrix::DataMatrix(double* input, size_t nrows, size_t ncols) : DataMatrix(nrows, ncols) {
  // copy data
  std::memcpy(this->data, input, nrows * ncols * sizeof(double));
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
      double value = std::atof(&(serializedVector[i]));
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
  double* newdata = new double[nrows * this->ncols];
  // copy entries of old matrix
  std::memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols * sizeof(double));
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
    double* newdata = new double[nrows * ncols];
    // copy entries of old matrix
    std::memcpy(newdata, this->data,
                std::min(this->nrows * this->ncols, nrows * ncols) * sizeof(double));
    delete[] this->data;
    this->data = newdata;
  }

  this->nrows = nrows;
  this->ncols = ncols;
  this->unused = 0;
}

void DataMatrix::resizeQuadratic(size_t size) {
  // Throw exception if matrix is not quadratic
  if (this->nrows != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::resizeQuadratic : DataMatrix is not quadratic");
  } else if (this->nrows == size) {
    // don't do anything, if matrix already has the correct size
    return;
  }

  // create new matrix
  double* newdata = new double[size * size];

  for (size_t i = 0; i < std::min(this->ncols, size); i++) {
    std::memcpy(&newdata[i * size], this->data + i * (this->ncols),
                std::min(this->ncols, size) * sizeof(double));
  }

  // Set new characteristics of DataMatrix
  delete[] this->data;
  this->data = newdata;
  this->nrows = size;
  this->ncols = size;
  this->unused = 0;
}

void DataMatrix::resizeZero(size_t nrows) {
  // don't do anything, if matrix already has the correct size
  if (nrows == this->nrows) {
    return;
  }

  // create new matrix
  double* newdata = new double[nrows * this->ncols];
  // copy entries of old matrix
  std::memcpy(newdata, this->data, std::min(this->nrows, nrows) * this->ncols * sizeof(double));

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
    double* newdata = new double[nrows * ncols];
    // copy entries of old matrix
    std::memcpy(newdata, this->data,
                std::min(this->nrows * this->ncols, nrows * ncols) * sizeof(double));

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

void DataMatrix::resizeToSubMatrix(size_t row_1, size_t col_1, size_t row_2, size_t col_2) {
  if ((row_1 > row_2) ||(col_1 > col_2)) {
    throw sgpp::base::data_exception(
        "DataMatrix::getSubMatrix : Expected indices do not fulfill requirements");
  } else if ((this->nrows < row_2) || (this->ncols < col_2)) {
    throw sgpp::base::data_exception("DataMatrix::getSubMatrix : Indices are out of bounds");
  } else if ((row_1 == row_2) && (col_1 == col_2)) {
    // do nothing
    return;
  }

  // create new matrix
  double* newdata = new double[(row_2 - row_1 + 1) * (col_2 - col_1 + 1)];

  for (size_t i = 0; i < (row_2 - row_1 + 1); i++) {
    std::memcpy(&newdata[i * (col_2 - col_1 + 1)],
                this->data + (row_1 - 1) * this->ncols + (col_1 -1) + i * this->ncols,
                (col_2 - col_1 + 1) * sizeof(double));
  }

  // Set new characteristics of DataMatrix
  delete[] this->data;
  this->data = newdata;
  this->nrows = row_2 - row_1 + 1;
  this->ncols = col_2 - col_1 + 1;
  this->unused = 0;
}

void DataMatrix::addSize(size_t inc_nrows) {
  // create new vector
  double* newdata = new double[(this->nrows + inc_nrows) * this->ncols];
  // copy entries of old vector
  std::memcpy(newdata, this->data, this->nrows * this->ncols * sizeof(double));

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
  double* newData = new double[nrows * ncols];

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
    throw sgpp::base::data_exception("DataMatrix::appendRow : Dimensions do not match");
  }

  size_t x = appendRow();
  // copy data
  std::memcpy(&this->data[x * this->ncols], vec.getPointer(), this->ncols * sizeof(double));
  return x;
}

size_t DataMatrix::appendCol(const DataVector& vec) {
  if (vec.getSize() != this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::appendCol : Dimensions do not match");
  }

  // create new vector
  double* newdata = new double[this->nrows  * (this->ncols + 1)];

  for (size_t i = 0; i < this->nrows; i++) {
    std::memcpy(&newdata[i * (this->ncols+1)], this->data + i * (this->ncols),
                this->ncols * sizeof(double));
  }

  for (size_t j = 0; j < this->nrows; j++) {
    newdata[j * (this->ncols + 1) + this->ncols] = vec[j];
  }

  delete[] this->data;
  this->data = newdata;
  this->ncols++;
  unused = 0;
  return(this->ncols);
}

void DataMatrix::setAll(double value) {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] = value;
  }
}

void DataMatrix::getRow(size_t row, DataVector& vec) const {
  if (vec.getSize() != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::getRow : Dimensions do not match");
  }

  for (size_t i = 0; i < this->ncols; i++) {
    vec[i] = this->data[row * ncols + i];
  }
}

void DataMatrix::getRow(size_t row, std::vector<double>& vec) const {
  vec.clear();

  for (size_t i = 0; i < this->ncols; i++) {
    vec.push_back(data[row * ncols + i]);
  }
}

void DataMatrix::setRow(size_t row, const DataVector& vec) {
  if (vec.getSize() != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::setRow : Dimensions do not match");
  }

  for (size_t i = 0; i < this->ncols; i++) {
    this->data[row * ncols + i] = vec.get(i);
  }
}

void DataMatrix::getColumn(size_t col, DataVector& vec) const {
  if (vec.getSize() != this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::getColumn : Dimensions do not match");
  }

  for (size_t j = 0; j < this->nrows; j++) {
    vec[j] = data[j * ncols + col];
  }
}

void DataMatrix::setColumn(size_t col, const DataVector& vec) {
  if (vec.getSize() != this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::setColumn : Dimensions do not match");
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
   this->data = new double[nrows * ncols];
   }
   */
  std::memcpy(this->data, matr.data,
              std::min(this->nrows * this->ncols, matr.nrows * matr.ncols) * sizeof(double));
}

/*
 void DataMatrix::copySmall(const DataMatrix& vec) {
 if (this == &vec) {
 return;
 }

 if (vec.ncols != 1 || ncols != 1 || nrows < vec.nrows) {
 return;
 }
 std::memcpy(this->data, vec.data, vec.nrows * sizeof(double));
 }

 DataMatrix& DataMatrix::operator=(const DataMatrix &vec) {
 if (this == &vec) {
 return *this;
 }

 if (nrows != vec.nrows || ncols != vec.ncols) {
 delete[] data;
 nrows = vec.nrows;
 ncols = vec.ncols;
 this->data = new double[nrows * ncols];
 }
 std::memcpy(this->data, vec.data, nrows * ncols * sizeof(double));
 return *this;
 }
 */
DataMatrix& DataMatrix::operator=(const DataMatrix& matr) {
  if (this == &matr) {
    return *this;
  }

  if (nrows * ncols != matr.ncols * matr.nrows) {
    // throw sgpp::base::data_exception("DataMatrix::= : Dimensions do not match");
    delete[] this->data;
    this->data = new double[matr.nrows * matr.ncols];
  }

  this->nrows = matr.nrows;
  this->ncols = matr.ncols;
  std::memcpy(this->data, matr.data, nrows * ncols * sizeof(double));
  return *this;
}

void DataMatrix::add(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::add : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] += matr.data[i];
  }
}

void DataMatrix::sub(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::sub : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] -= matr.data[i];
  }
}

void DataMatrix::addReduce(DataVector& reduction) {
  if (this->nrows != reduction.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::addReduce : Dimensions do not match");
  }

  for (size_t i = 0; i < this->nrows; i++) {
    double tmp = 0.0;

    for (size_t j = 0; j < this->ncols; j++) {
      tmp += this->data[(i * this->ncols) + j];
    }

    reduction.set(i, tmp);
  }
}

void DataMatrix::addReduce(DataVector& reduction, DataVector& beta, size_t start_beta) {
  if (this->nrows != reduction.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::addReduce : Dimensions do not match (reduction)");
  }

  if (this->ncols + start_beta > beta.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::addReduce : Dimensions do not match (beta)");
  }

  for (size_t i = 0; i < this->nrows; i++) {
    double tmp = 0.0;

    for (size_t j = 0; j < this->ncols; j++) {
      tmp += beta[j + start_beta] * this->data[(i * this->ncols) + j];
    }

    reduction.set(i, reduction[i] + tmp);
  }
}

void DataMatrix::expand(const DataVector& expand) {
  if (this->nrows != expand.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::expand : Dimensions do not match");
  }

  for (size_t i = 0; i < this->nrows; i++) {
    for (size_t j = 0; j < this->ncols; j++) {
      this->data[(i * this->ncols) + j] = expand.get(i);
    }
  }
}

void DataMatrix::componentwise_mult(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::componentwise_mult : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] *= matr.data[i];
  }
}

void DataMatrix::componentwise_div(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::componentwise_div : Dimensions do not match");
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

 void DataMatrix::getLine(int row, std::vector<double>& vec) {
 vec.clear();

 for (int i = 0; i < this->ncols; i++) {
 vec.push_back(data[row * ncols + i]);
 }
 }
 */

/*
 double DataMatrix::dotProduct(DataMatrix &vec) {
 double sum = 0.0;

 for (int i = 0; i < nrows; i++) {
 sum += data[i] * vec.data[i];
 }
 return sum;
 }
 */

void DataMatrix::mult(double scalar) {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; i++) {
    data[i] *= scalar;
  }
}

void DataMatrix::mult(const DataVector& x, DataVector& y) {
  if (ncols != x.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::mult : Dimensions do not match (x)");
  }

  if (nrows != y.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::mult : Dimensions do not match (y)");
  }

  for (size_t i = 0; i < nrows; i++) {
    double entry = 0.0;

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

double DataMatrix::sum() const {
  size_t n = nrows * ncols;
  double result = 0.0;

  for (size_t i = 0; i < n; i++) {
    result += data[i];
  }

  return result;
}
/*
 double DataMatrix::maxNorm() {
 int n = nrows * ncols;
 double max = 0.0;
 for (int i = 0; i < n; i++) {
 if (max < fabs(data[i]))
 {
 max = fabs(data[i]);
 }
 }
 return max;
 }

 void DataMatrix::partitionClasses(double border) {
 int n = nrows * ncols;
 for (int i = 0; i < n; i++) {
 data[i] = data[i] > border ? 1.0 : -1.0;
 }
 }

 void DataMatrix::axpy(double alpha, DataMatrix& x) {
 if (nrows != x.nrows || ncols != x.ncols) {
 return;
 }
 int n = nrows * ncols;
 double* p_x = x.data;
 double* p_d = data;

 for (int i = 0; i < n; i++) {
 p_d[i] += alpha * p_x[i];
 }
 }
 */

void DataMatrix::normalizeDimension(size_t d) { normalizeDimension(d, 0.0); }

void DataMatrix::normalizeDimension(size_t d, double border) {
  size_t n = nrows * ncols;

  if (ncols <= d) {
    throw sgpp::base::data_exception(
        "DataMatrix::normalizeDimension : Not enough columns in DataMatrix");
  }

  // determine min and max
  double xmin, xmax;
  minmax(d, &xmin, &xmax);

  double delta = (xmax - xmin) / (1 - 2 * border);

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

double DataMatrix::min(size_t d) const {
  size_t n = nrows * ncols;
  double min = INFINITY;

  for (size_t i = d; i < n; i += ncols) {
    if (min > data[i]) {
      min = data[i];
    }
  }

  return min;
}

double DataMatrix::min() const {
  size_t n = nrows * ncols;
  double min = INFINITY;

  for (size_t i = 0; i < n; i++) {
    if (min > data[i]) {
      min = data[i];
    }
  }

  return min;
}

double DataMatrix::max(size_t d) const {
  size_t n = nrows * ncols;
  double max = -INFINITY;

  for (size_t i = d; i < n; i += ncols) {
    if (max < data[i]) {
      max = data[i];
    }
  }

  return max;
}

double DataMatrix::max() const {
  size_t n = nrows * ncols;
  double max = -INFINITY;

  for (size_t i = 0; i < n; i++) {
    if (max < data[i]) {
      max = data[i];
    }
  }

  return max;
}

void DataMatrix::minmax(size_t col, double* min, double* max) const {
  size_t n = nrows * ncols;

  if (ncols <= col) {
    throw sgpp::base::data_exception("DataMatrix::minmax : Not enough entries in DataMatrix");
  }

  // find min and max of column col
  double min_t = INFINITY;
  double max_t = -INFINITY;

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

void DataMatrix::minmax(double* min, double* max) const {
  size_t n = nrows * ncols;

  double min_t = INFINITY;
  double max_t = -INFINITY;

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

double* DataMatrix::getPointer() { return data; }

const double* DataMatrix::getPointer() const { return data; }

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
}  // namespace sgpp
