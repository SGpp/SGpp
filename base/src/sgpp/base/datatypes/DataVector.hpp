// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAVECTOR_HPP
#define DATAVECTOR_HPP

#include <sgpp/globaldef.hpp>

#include <initializer_list>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * A class to store one-dimensional data.
 * Typically, an object of type DataVector will contain an array
 * of (hierarchical) coefficients (or surplusses), or the coordinates
 * of a data point at which a sparse grid function should be
 * evaluated.
 */
class DataVector : public std::vector<double> {
 public:
  /**
   * Create an empty DataVector.
   */
  DataVector();

  /**
   * Copy constructor.
   */
  DataVector(const DataVector&) = default;

  /**
   * Move constructor
   */
  DataVector(DataVector&&) = default;

  /**
   * Copy assignment operator
   */
  DataVector& operator=(const DataVector&) = default;

  /**
   * Move assignment operator
   */
  DataVector& operator=(DataVector&&) = default;

  /**
   * Destructor
   */
  ~DataVector() = default;

  /**
   * Create a DataVector with @em size elements (uninitialized values).
   *
   * @param size Number of elements
   */
  explicit DataVector(size_t size);

  /**
   * Create a DataVector with @em size elements and initializes
   * all elements with the same value.
   *
   * @param size Number of elements
   * @param value Value for all entries
   */
  DataVector(size_t size, double value);

  /**
   * Create a new DataVector from a double array with size elements.
   *
   * @param input double array that contains the data
   * @param size number of elements
   */
  DataVector(double* input, size_t size);

  /**
   * Create a new DataVector from a std::vector<double>.
   *
   * @param input std::vector<double> that contains the data
   */
  explicit DataVector(std::vector<double> input);

  /**
   * Create a new DataVector from a std::initializer_list<double>.
   *
   * @param input std::initializer_list<double> that contains the data
   */
  explicit DataVector(std::initializer_list<double> input);

  /**
   * Create a new DataVector from a std::vector<int>.
   *
   * @param input std::vector<int> that contains the data
   */
  explicit DataVector(std::vector<int> input);

  static DataVector fromFile(const std::string& fileName);

  static DataVector fromString(const std::string& serializedVector);

  /**
   * Resizes the DataVector to size elements.
   * All new additional entries are set to zero.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   *
   * @param nrows New number of rows of the DataVector
   */
  void resizeZero(size_t nrows);

  /**
   * Resizes the DataVector by removing entries. Throws an exception
   * if boundaries a violated.
   *
   * @param remainingIndex vector that contains the remaining indices of the DataVector
   */
  void restructure(std::vector<size_t>& remainingIndex);

  /**
   * Removes indexes form the vector. Throws an exception if the boundaries are violated
   *
   * @param indexesToRemove a vector if indexes that will be removed from the vector
   */
  void remove(std::vector<size_t>& indexesToRemove);

  /**
   * Appends a new element and returns index of it.
   *
   * @return Index of new element
   */
  size_t append();

  /**
   * Appends a new element and returns index of new element.
   *
   * @param value Value of new element
   * @return Index of new element
   */
  size_t append(double value);

  /**
   * Appends sequence of elements from another DataVector given by iterators
   * @param first	iterator pointing to the first element of the sequence
   * @param last 	iterator pointing to the last element of the sequence
   */
  void append(DataVector::iterator first, DataVector::iterator last);

  /**
   * Inserts a new element at the given index.
   *
   * @param index Index of new element
   * @param value Value of new element
   */
  void insert(size_t index, double value);

  /**
   * Sets all values of DataVector to value
   *
   * @param value New value for all entries
   */
  void setAll(double value);

  /**
   * Returns the i-th element.
   *
   * @param i position of the element
   * @return data[i]
   */
  inline double get(size_t i) const { return (*this)[i]; }

  /**
   * Sets the element at index i to value.
   *
   * @param i Index
   * @param value New value for element
   */
  void set(size_t i, double value);

  /**
   * Copies the data from another DataVector vec.
   * Disregards the number of entries set for the two vectors,
   * i.e., just copies the data entry by entry.
   * If the size matches, the current DataVector is an
   * exact copy of vec. If not, as many elements as possible are
   * copied, and everything else is left untouched.
   *
   * @param vec The source DataVector containing the data
   */
  void copyFrom(const DataVector& vec);

  /**
   * Adds the values from another DataVector to the current values.
   * Modifies the current values.
   *
   * @param vec The DataVector which is added to the current values
   */
  void add(const DataVector& vec);

  /***
   * Accumulation (summation) of vectors using Kahan's summation formula
   * for better precision.
   *
   * All vectors need to be added to one for the summation formula to work. It
   * will not work in a tree-like summation pattern.
   *
   * @param vec The DataVector that will be added to the current one
   */
  void accumulate(const DataVector& vec);

  /**
   * Subtracts the values from another DataVector of the current values.
   * Modifies the current values.
   *
   * @param vec The DataVector which is subtracted from the current values
   */
  void sub(const DataVector& vec);

  /**
   * Multiplies the current DataVector component-wise with another DataVector.
   * Modifies the current values.
   * Performs
   * @code
   * for i from 1 to this.getSize()
   *   this[i] *= vec[i]
   * @endcode
   *
   * @param vec the DataVector which is multiplied to current DataVector
   */
  void componentwise_mult(const DataVector& vec);

  /**
   * Divides the current DataVector component-wise by another DataVector.
   * Modifies the current values.
   * Performs
   * @code
   * for i from 1 to this.getSize()
   *   this[i] /= vec[i]
   * @endcode
   * Note: <b>No check for division by zero!</b>
   *
   * @param vec the DataVector which the current DataVector is divided by
   */
  void componentwise_div(const DataVector& vec);

  /**
   * Returns the dot product of the two vectors.
   *
   * @param vec Reference to another vector
   *
   * @return The dot-product
   */
  double dotProduct(const DataVector& vec) const;

  /**
   * multiplies all elements by a constant factor
   *
   * @param scalar the constant
   */
  void mult(double scalar);

  /**
   * Squares all elements of the DataVector
   */
  void sqr();

  /**
   * Takes the square root of all elements of the DataVector
   */
  void sqrt();

  /**
   * Sets all elements to their absolute value.
   *
   */
  void abs();

  /**
   * Returns the sum of all elements
   *
   * @return The sum of all elements
   */
  double sum() const;

  /**
   * calculates the vector's max norm
   *
   * @return the vector's max norm
   */
  double maxNorm() const;

  /**
   * Returns the vector's root mean square (RMS)-norm, i.e.,
   * @f$\sqrt{ 1/N \sum_{i=1}^N x_i^2 }@f$. If the vector's entries
   * correspond to function values on a full grid, this is the
   * discrete @f$L^2@f$-norm of the corresponding function.
   *
   * @return The vector's root mean square-norm.
   */
  double RMSNorm() const;

  /**
   * Returns the vector's @f$l^2@f$-norm, i.e.,
   * @f$\sqrt{ \sum_i x_i^2 }@f$.
   *
   * @return The vector's @f$l^2@f$-norm.
   */
  double l2Norm() const;

  /**
   * Returns the minimum over all entries.
   *
   * @return Minimal value
   */
  double min() const;

  /**
   * Returns the maximum over all entries.
   *
   * @return global maximum
   */
  double max() const;

  /**
   * Determines minimum and maximum over all entries.
   *
   * @param min Reference variable for the minimum
   * @param max Reference variable for the maximum
   */
  void minmax(double* min, double* max) const;

  /**
   * Adds a*x to current vector.
   * BLAS Level 1 (elementary vector operations) operation: axpy.
   *
   * @param a A scalar
   * @param x Reference to the DataVector
   */
  void axpy(double a, DataVector& x);

  /**
   * gets a pointer to the data array
   *
   * @return pointer to the data array
   */
  double* getPointer();

  /**
   * gets a const pointer to the data array
   *
   * @return const pointer to the data array
   */
  const double* getPointer() const;

  /**
   * gets the elements stored in the vector
   * \deprecated in favour of the equivalent size() method
   *
   * @return elements stored in the vector
   */
  inline size_t getSize() const { return this->size(); }

  /**
   * Determines the number of non-zero elements in the vector.
   *
   * @return The number of non-zero elements
   */
  size_t getNumberNonZero() const;

  /**
   * Partitions vector into two classes using a choosen border.
   *
   * @param threshold value of the border
   */
  void partitionClasses(double threshold);

  /**
   * Normalizes vector entries to [0,1]
   *
   */
  void normalize();

  /**
   * Normalizes vector entries to [border, 1-border]
   *
   * @param border width of border
   */
  void normalize(double border);

  /**
   * Writes the data stored in the DataVector into a string
   *
   * @param text string to which the data is written
   */
  void toString(std::string& text) const;

  /**
   * Returns a description of the DataVector as a string.
   *
   * @returns string of the DataVector
   */
  std::string toString() const;

  void toFile(const std::string& fileName) const;

 private:
  using std::vector<double>::insert;
  /// Corrections for Kahan's summation in accumulate()
  std::vector<double> correction;
};

}  // namespace base
}  // namespace sgpp

#endif /* DATAVECTOR_HPP */
