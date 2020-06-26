// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

DataMatrix::DataMatrix() : DataMatrix(0, 0) {}

DataMatrix::DataMatrix(size_t nrows, size_t ncols) : DataMatrix(nrows, ncols, 0.0) {}

DataMatrix::DataMatrix(size_t nrows, size_t ncols, double value) : nrows(nrows), ncols(ncols) {
  this->assign(nrows * ncols, value);
}

DataMatrix::DataMatrix(const double* input, size_t nrows, size_t ncols)
    : std::vector<double>(input, input + nrows * ncols), nrows(nrows), ncols(ncols) {}

DataMatrix::DataMatrix(std::vector<double> input, size_t nrows)
    : DataMatrix(input.data(), nrows, input.size() / nrows) {}

DataMatrix::DataMatrix(std::initializer_list<double> input, size_t nrows)
    : std::vector<double>(input), nrows(nrows), ncols(input.size() / nrows) {}

DataMatrix DataMatrix::fromFile(const std::string& fileName) {
  std::ifstream f(fileName, std::ifstream::in);
  f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  std::string content;
  content.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
  return DataMatrix::fromString(content);
}

DataMatrix DataMatrix::fromString(const std::string& serializedVector) {
  DataMatrix m(0, 0);

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
      ++i;
    } else if (state == PARSER_STATE::ROWVALUE) {
      //      size_t next;
      double value = std::atof(&(serializedVector[i]));
      row.append(value);
      state = PARSER_STATE::ROWCOMMAEND;
      //      i += next;
      while (serializedVector[i] != ',' && serializedVector[i] != ']') ++i;
    } else if (state == PARSER_STATE::ROWCOMMAEND) {
      if (c == ',') {
        state = PARSER_STATE::ROWVALUE;
        ++i;
      } else if (c == ']') {
        if (m.getNcols() == 0 || m.getNrows() == 0) {
          // set up the dimension after having read the first row
          m.resize(0, row.getSize());
        }
        m.appendRow(row);
        state = PARSER_STATE::COMMAEND;
        ++i;
      }
    } else if (state == PARSER_STATE::COMMAEND) {
      if (c == ',') {
        state = PARSER_STATE::ROW;
        ++i;
      } else if (c == ']') {
        state = PARSER_STATE::END;
        ++i;
      }
    } else if (state == PARSER_STATE::END) {
      // only reached if a non-whitespace character was encountered after closing brace
      throw data_exception("error: could not parse DataMatrix file");
    }
  }
  return m;
}

void DataMatrix::resize(size_t nrows) { this->resizeRows(nrows); }

void DataMatrix::resizeRows(size_t nrows) {
  // don't do anything, if matrix already has the correct size
  if (nrows == this->nrows) {
    return;
  }
  this->nrows = nrows;
  this->std::vector<double>::resize(nrows * ncols);
}

void DataMatrix::resize(size_t nrows, size_t ncols) { this->resizeRowsCols(nrows, ncols); }

void DataMatrix::resizeRowsCols(size_t nrows, size_t ncols) {
  // don't do anything, if matrix already has the correct size
  if ((nrows == this->nrows) && (ncols == this->ncols)) {
    return;
  }
  this->nrows = nrows;
  this->ncols = ncols;
  this->std::vector<double>::resize(nrows * ncols);
}

void DataMatrix::resizeQuadratic(size_t size) {
  // Throw exception if matrix is not quadratic
  if (this->nrows != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::resizeQuadratic : DataMatrix is not quadratic");
  } else if (this->nrows == size) {
    // don't do anything, if matrix already has the correct size
    return;
  }

  DataMatrix oldMatrix(*this);

  size_t min_size = std::min(this->ncols, size);
  this->resize(size, size);
  for (size_t i = 0; i < min_size; ++i) {
    std::copy(oldMatrix.row_begin(i), oldMatrix.row_begin(i) + min_size, this->row_begin(i));
  }
}

void DataMatrix::resizeZero(size_t nrows) { this->resizeRows(nrows); }

void DataMatrix::resizeZero(size_t nrows, size_t ncols) { this->resizeRowsCols(nrows, ncols); }

void DataMatrix::resizeToSubMatrix(size_t row_1, size_t col_1, size_t row_2, size_t col_2) {
  if ((row_1 > row_2) || (col_1 > col_2)) {
    throw sgpp::base::data_exception(
        "DataMatrix::getSubMatrix : Expected indices do not fulfill requirements");
  } else if ((this->nrows < row_2) || (this->ncols < col_2)) {
    throw sgpp::base::data_exception("DataMatrix::getSubMatrix : Indices are out of bounds");
  } else if ((row_1 == row_2) && (col_1 == col_2)) {
    // do nothing
    return;
  }

  auto regionBegin = this->begin() + row_1 * this->ncols + col_1;
  auto regionEnd = this->begin() + row_2 * this->ncols + col_2;
  size_t ncols_new = col_2 - col_1 + 1;

  // copy region row-wise
  size_t resizedOffset = 0;
  for (auto it = regionBegin; it < regionEnd + 1; it += this->ncols) {
    std::copy(it, it + (col_2 - col_1) + 1, this->begin() + resizedOffset);
    resizedOffset += ncols_new;
  }

  this->nrows = row_2 - row_1 + 1;
  this->ncols = ncols_new;

  // free unused memory
  this->std::vector<double>::resize(this->nrows * this->ncols);
  this->shrink_to_fit();
}

void DataMatrix::reserveAdditionalRows(size_t inc_nrows) {
  this->reserve(this->size() + inc_nrows * this->ncols);
}

size_t DataMatrix::appendRow() {
  this->insert(this->end(), this->ncols, 0.0);

  size_t x = nrows;
  ++this->nrows;

  return x;
}

void DataMatrix::transpose() {
  if (this->nrows == this->ncols) {
    for (size_t i = 1; i < this->nrows; ++i) {
      for (size_t j = 0; j < i; ++j) {
        std::swap((*this)[(j * this->nrows) + i], (*this)[(i * this->ncols) + j]);
      }
    }
  } else {
    DataMatrix oldMatrix(*this);
    this->resize(this->ncols, this->nrows);
    for (size_t i = 0; i < oldMatrix.nrows; ++i) {
      for (size_t j = 0; j < oldMatrix.ncols; ++j) {
        (*this)(j, i) = oldMatrix(i, j);
      }
    }
  }
}

size_t DataMatrix::appendRow(const DataVector& vec) {
  if (vec.getSize() != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::appendRow : Dimensions do not match");
  }
  this->insert(this->end(), vec.begin(), vec.end());
  ++this->nrows;
  return this->nrows - 1;
}

size_t DataMatrix::appendCol(const DataVector& vec) {
  if (vec.getSize() != this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::appendCol : Dimensions do not match");
  }
  if (this->nrows == 0) {
    ++this->ncols;
    return this->ncols - 1;
  }
  if (this->nrows == 1) {
    this->push_back(vec[0]);
    ++this->ncols;
    return this->ncols - 1;
  }
  if (this->ncols == 0) {
    this->assign(vec.begin(), vec.end());
    ++this->ncols;
    return this->ncols - 1;
  }

  this->reserve(this->size() + this->nrows);
  size_t initial_size = this->size();
  size_t new_ncols = this->ncols + 1;

  size_t retained_rows = initial_size / new_ncols;
  size_t retained_rows_rest = initial_size % new_ncols;

  // append partial row
  this->insert(this->end(), this->row_begin(retained_rows) + retained_rows_rest,
               this->row_end(retained_rows));
  this->push_back(vec[retained_rows]);

  // append full rows
  for (size_t push_back_row = retained_rows + 1; push_back_row < this->nrows; ++push_back_row) {
    this->insert(this->end(), this->row_begin(push_back_row), this->row_end(push_back_row));
    this->push_back(vec[push_back_row]);
  }

  // copy partial row (if anything is left to do)
  std::copy_backward(this->row_begin(retained_rows),
                     this->row_begin(retained_rows) + retained_rows_rest,
                     this->begin() + retained_rows * new_ncols + retained_rows_rest);
  (*this)[retained_rows * new_ncols - 1] = vec[retained_rows - 1];

  // copy rest of the rows
  for (size_t copy_row = retained_rows - 1; copy_row > 0; --copy_row) {
    std::copy_backward(this->row_begin(copy_row), this->row_end(copy_row),
                       this->begin() + (copy_row + 1) * new_ncols - 1);
    (*this)[copy_row * new_ncols - 1] = vec[copy_row - 1];
  }

  this->ncols = new_ncols;
  return this->ncols - 1;
}

void DataMatrix::setAll(double value) {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] = value;
  }
}

void DataMatrix::getRow(size_t row, DataVector& vec) const {
  if (vec.getSize() != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::getRow : Dimensions do not match");
  }

  for (size_t i = 0; i < this->ncols; ++i) {
    vec[i] = (*this)[row * ncols + i];
  }
}

void DataMatrix::getRow(size_t row, std::vector<double>& vec) const {
  vec.clear();

  for (size_t i = 0; i < this->ncols; ++i) {
    vec.push_back((*this)[row * ncols + i]);
  }
}

void DataMatrix::setRow(size_t row, const DataVector& vec) {
  if (vec.getSize() != this->ncols) {
    throw sgpp::base::data_exception("DataMatrix::setRow : Dimensions do not match");
  } else if (row >= this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::setRow : \"row\" out of bounds");
  }

  for (size_t i = 0; i < this->ncols; ++i) {
    (*this)[row * ncols + i] = vec.get(i);
  }
}

void DataMatrix::getColumn(size_t col, DataVector& vec) const {
  if (vec.getSize() != this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::getColumn : Dimensions do not match");
  }

  for (size_t j = 0; j < this->nrows; ++j) {
    vec[j] = (*this)[j * ncols + col];
  }
}

void DataMatrix::setColumn(size_t col, const DataVector& vec) {
  if (vec.getSize() != this->nrows) {
    throw sgpp::base::data_exception("DataMatrix::setColumn : Dimensions do not match");
  }

  for (size_t j = 0; j < this->nrows; ++j) {
    (*this)[j * ncols + col] = vec.get(j);
  }
}

void DataMatrix::copyFrom(const DataMatrix& matr) {
  if (*this == matr) {
    return;
  }
  std::copy(matr.begin(), matr.begin() + std::min(this->size(), matr.size()), this->begin());
}

void DataMatrix::add(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::add : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] += matr[i];
  }
}

void DataMatrix::sub(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::sub : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] -= matr[i];
  }
}

void DataMatrix::addReduce(DataVector& reduction) {
  if (this->nrows != reduction.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::addReduce : Dimensions do not match");
  }

  for (size_t i = 0; i < this->nrows; ++i) {
    double tmp = 0.0;

    for (size_t j = 0; j < this->ncols; ++j) {
      tmp += (*this)[(i * this->ncols) + j];
    }

    reduction.set(i, tmp);
  }
}

void DataMatrix::addReduce(DataVector& reduction, DataVector& beta, size_t start_beta) {
  if (this->nrows != reduction.getSize()) {
    throw sgpp::base::data_exception(
      "DataMatrix::addReduce : Dimensions do not match (reduction)");
  }

  if (this->ncols + start_beta > beta.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::addReduce : Dimensions do not match (beta)");
  }

  for (size_t i = 0; i < this->nrows; ++i) {
    double tmp = 0.0;

    for (size_t j = 0; j < this->ncols; ++j) {
      tmp += beta[j + start_beta] * (*this)[(i * this->ncols) + j];
    }

    reduction.set(i, reduction[i] + tmp);
  }
}

void DataMatrix::expand(const DataVector& expand) {
  if (this->nrows != expand.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::expand : Dimensions do not match");
  }

  for (size_t i = 0; i < this->nrows; ++i) {
    for (size_t j = 0; j < this->ncols; ++j) {
      (*this)[(i * this->ncols) + j] = expand.get(i);
    }
  }
}

void DataMatrix::componentwise_mult(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::componentwise_mult : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] *= matr[i];
  }
}

void DataMatrix::componentwise_div(const DataMatrix& matr) {
  if (this->nrows != matr.nrows || this->ncols != matr.ncols) {
    throw sgpp::base::data_exception("DataMatrix::componentwise_div : Dimensions do not match");
  }

  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] /= matr[i];
  }
}

void DataMatrix::mult(double scalar) {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] *= scalar;
  }
}

void DataMatrix::mult(const DataVector& x, DataVector& y) const {
  if (ncols != x.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::mult : Dimensions do not match (x)");
  }

  if (nrows != y.getSize()) {
    throw sgpp::base::data_exception("DataMatrix::mult : Dimensions do not match (y)");
  }

  for (size_t i = 0; i < nrows; ++i) {
    double entry = 0.0;

    for (size_t j = 0; j < ncols; ++j) {
      entry += (*this)[(i * ncols) + j] * x[j];
    }

    y[i] = entry;
  }
}

void DataMatrix::sqr() {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] = (*this)[i] * (*this)[i];
  }
}

void DataMatrix::sqrt() {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] = std::sqrt((*this)[i]);
  }
}

void DataMatrix::abs() {
  size_t n = nrows * ncols;

  for (size_t i = 0; i < n; ++i) {
    (*this)[i] = std::abs((*this)[i]);
  }
}

double DataMatrix::sum() const {
  size_t n = nrows * ncols;
  double result = 0.0;

  for (size_t i = 0; i < n; ++i) {
    result += (*this)[i];
  }

  return result;
}

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

  double delta = (xmax - xmin) / (1.0 - 2.0 * border);

  if (delta == 0.0) {
    for (size_t i = d; i < n; i += ncols) {
      (*this)[i] = 0.5;
    }
  } else {
    for (size_t i = d; i < n; i += ncols) {
      (*this)[i] = ((*this)[i] - xmin) / delta + border;
    }
  }
}

void DataMatrix::toString(std::string& text) const {
  std::stringstream str;

  str << std::scientific;
  str.precision(20);

  str << "[";

  for (size_t i = 0; i < nrows; ++i) {
    str << "[";

    for (size_t j = 0; j < ncols; ++j) {
      if (j != 0) {
        str << ", ";
        // add linebreak for readability
        if (j % 20 == 0) {
          str << std::endl;
        }
      }

      str << (*this)[i * ncols + j];
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
  double min = std::numeric_limits<double>::infinity();

  for (size_t i = d; i < n; i += ncols) {
    if (min > (*this)[i]) {
      min = (*this)[i];
    }
  }

  return min;
}

double DataMatrix::min() const {
  size_t n = nrows * ncols;
  double min = std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < n; ++i) {
    if (min > (*this)[i]) {
      min = (*this)[i];
    }
  }

  return min;
}

double DataMatrix::max(size_t d) const {
  size_t n = nrows * ncols;
  double max = -std::numeric_limits<double>::infinity();

  for (size_t i = d; i < n; i += ncols) {
    if (max < (*this)[i]) {
      max = (*this)[i];
    }
  }

  return max;
}

double DataMatrix::max() const {
  size_t n = nrows * ncols;
  double max = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < n; ++i) {
    if (max < (*this)[i]) {
      max = (*this)[i];
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
  double min_t = std::numeric_limits<double>::infinity();
  double max_t = -std::numeric_limits<double>::infinity();

  for (size_t i = col; i < n; i += ncols) {
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

void DataMatrix::minmax(double* min, double* max) const {
  size_t n = nrows * ncols;

  double min_t = std::numeric_limits<double>::infinity();
  double max_t = -std::numeric_limits<double>::infinity();

  for (size_t i = 0; i < n; ++i) {
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

double* DataMatrix::getPointer() { return const_cast<double*>(this->data()); }

const double* DataMatrix::getPointer() const { return this->data(); }

size_t DataMatrix::getNumberNonZero() const {
  size_t n = nrows * ncols;
  size_t nonZero = 0;

  for (size_t i = 0; i < n; ++i) {
    if (std::abs((*this)[i]) > 0.0) {
      ++nonZero;
    }
  }

  return nonZero;
}

}  // namespace base
}  // namespace sgpp
