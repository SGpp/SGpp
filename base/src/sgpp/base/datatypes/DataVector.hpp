// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAVECTOR_HPP
#define DATAVECTOR_HPP

#include <string>
#include <vector>

#include <sgpp/base/datatypes/DataVectorDefinition.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * A class to store one-dimensional data.
     * Typically, an object of type DataVector will contain an array
     * of (hierarchical) coefficients (or surplusses), or the coordinates
     * of a data point at which a sparse grid function should be
     * evaluated.
     */
    class DataVector {
      public:
        /**
         * Create a DataVector with @em size elements (uninitialized values).
         *
         * @param size Number of elements
         */
        DataVector(size_t size);

        /**
         * Create a DataVector with @em size elements and initializes
         * all elements with the same value.
         *
         * @param size Number of elements
         * @param value Value for all entries
         */
        DataVector(size_t size, float_t value);

        /**
         * Create a new DataVector that is a copy of vec.
         *
         * @param vec Reference to another instance of DataVector
         */
        DataVector(const DataVector& vec);

        /**
         * Create a new DataVector from a float_t array with size elements.
         *
         * @param input float_t array that contains the data
         * @param size number of elements
         */
        DataVector(float_t* input, size_t size);

        /**
         * Create a new DataVector from a std::vector<float_t>.
         *
         * @param input std::vector<float_t> that contains the data
         */
        DataVector(std::vector<float_t> input);

        /**
         * Create a new DataVector from a std::vector<int>.
         *
         * @param input std::vector<int> that contains the data
         */
        DataVector(std::vector<int> input);

        /**
         * Constructor that constructs a DataVector from a DataVectorDefinition structure
         *
         * @param DataVectorDef reference to a DataVectorDefinition structure
         */
        DataVector(DataVectorDefinition& DataVectorDef);

        /**
         * Fill DataVectorDefintion with data from the current DataVector
         *
         * @param DataVectorDef reference to a DataVectorDefinition struct
         */
        void getDataVectorDefinition(DataVectorDefinition& DataVectorDef);
        /**
         * Set DataVectorDefintion for current DataVector
         *
         * @param DataVectorDef reference to a DataVectorDefinition struct
         */
        void setDataVectorDefinition(DataVectorDefinition& DataVectorDef);

        /**
         * Resizes the DataVector to size elements.
         * All new additional entries are uninitialized.
         * If nrows is smaller than the current number of rows,
         * all superfluous entries are removed.
         *
         * @param size New number of elements of the DataVector
         */
        void resize(size_t size);

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
         * Add add potentially new elements to the DataVector. The size remains unchanged
         * Reserves memory for potentially inc_elems new elements;
         * the actual number of elements remains unchanged.
         * Corresponds to a resize to size+inc_elems new elements while leaving
         * the current vector's size unchanged.
         *
         * @param inc_elems Number of additional elements for which storage is to be reserved.
         */
        void addSize(size_t inc_elems);

        /**
         * Appends a new element and returns index of it.
         * If the new element does not fit into the reserved memory,
         * reserves memory for getInc() additional elements.
         *
         * @return Index of new element
         */
        size_t append();

        /**
         * Appends a new element and returns index of new element.
         * If the new element does not fit into the reserved memory,
         * reserves memory for getInc() additional elements.
         *
         * @param value Value of new element
         * @return Index of new element
         */
        size_t append(float_t value);

        /**
         * Inserts a new element at the given index.
         * If the new element does not fit into the reserved memory,
         * reserves memory for getInc() additional elements.
         *
         * @param index Index of new element
         * @param value Value of new element
         */
        void insert(size_t index, float_t value);

        /**
         * Sets all values of DataVector to value
         *
         * @param value New value for all entries
         */
        void setAll(float_t value);

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
         * Copies the data from another, smaller DataVector vec.
         * Has no effect, if vec has more elements than the current
         * vector. If it has the same or smaller size, the first
         * vec.size() elements of the current vector are overwritten.
         *
         * @param vec the DataVector containing the data
         * @return *this
         */
        //  void copySmall(const DataVector& vec);
        /**
         * Copies the data from another DataVector.
         * Dimensions have to match.
         *
         * @param vec the DataVector containing the data
         */
        DataVector& operator=(const DataVector& vec);

        /**
         * Returns a reference to the i-th element.
         *
         * @param i position of the element
         * @return data[i]
         */
        inline float_t& operator[](size_t i) {
          return data[i];
        }

        /**
         * Returns a constant reference to the i-th element.
         *
         * @param i position of the element
         * @return data[i]
         */
        inline const float_t& operator[](size_t i) const {
          return data[i];
        }

        /**
         * Returns the i-th element.
         *
         * @param i position of the element
         * @return data[i]
         */
        inline float_t get(size_t i) const {
          return data[i];
        }

        /**
         * Sets the element at index i to value.
         *
         * @param i Index
         * @param value New value for element
         */
        void set(size_t i, float_t value);

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
         * multiplies all elements by a constant factor
         *
         * @param scalar the constant
         */
        void mult(float_t scalar);

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
        float_t sum() const;

        /**
         * calculates the vector's max norm
         *
         * @return the vector's max norm
         */
        float_t maxNorm() const;

        /**
         * Returns the vector's root mean square (RMS)-norm, i.e.,
         * @f$\sqrt{ 1/N \sum_{i=1}^N x_i^2 }@f$. If the vector's entries
         * correspond to function values on a full grid, this is the
         * discrete @f$L^2@f$-norm of the corresponding function.
         *
         * @return The vector's root mean square-norm.
         */
        float_t RMSNorm() const;

        /**
         * Returns the vector's @f$l^2@f$-norm, i.e.,
         * @f$\sqrt{ \sum_i x_i^2 }@f$.
         *
         * @return The vector's @f$l^2@f$-norm.
         */
        float_t l2Norm() const;

        /**
         * Returns the minimum over all entries.
         *
         * @return Minimal value
         */
        float_t min() const;

        /**
         * Returns the maximum over all entries.
         *
         * @return global maximum
         */
        float_t max() const;

        /**
         * Determines minimum and maximum over all entries.
         *
         * @param min Reference variable for the minimum
         * @param max Reference variable for the maximum
         */
        void minmax(float_t* min, float_t* max) const;

        /**
         * Adds a*x to current vector.
         * BLAS Level 1 (elementary vector operations) operation: axpy.
         *
         * @param a A scalar
         * @param x Reference to the DataVector
         */
        void axpy(float_t a, DataVector& x);

        /**
         * Returns the dot product of the two vectors.
         *
         * @param vec Reference to another vector
         *
         * @return The dot-product
         */
        float_t dotProduct(const DataVector& vec) const;

        /**
         * gets a pointer to the data array
         *
         * @return pointer to the data array
         */
        float_t* getPointer();

        /**
         * gets a const pointer to the data array
         *
         * @return const pointer to the data array
         */
        const float_t* getPointer() const;

        /**
         * gets the elements stored in the vector
         *
         * @return elements stored in the vector
         */
        inline size_t getSize() const {
          return size;
        }

        /**
         * Returns the number of unused elements.
         *
         * @return number of unused elements
         */
        inline size_t getUnused() const {
          return unused;
        }

        /**
         * Determines the number of non-zero elements in the vector.
         *
         * @return The number of non-zero elements
         */
        size_t getNumberNonZero() const;

        /**
         * Get the current number of elements by which the DataVector is extended,
         * if append() or insert() is called and no unused rows are left
         *
         * @return Increment
         */
        inline size_t getInc() const {
          return inc_elems;
        }

        /**
         * Sets the current number of elements by which the DataVector is extended,
         * if append() or insert() is called and no unused elements are left.
         * Defaults to 100.
         *
         * @param inc_elems Increment
         */
        void setInc(size_t inc_elems) {
          this->inc_elems = inc_elems;
        }

        /**
         * Partitions vector into two classes using a choosen border.
         *
         * @param threshold value of the border
         */
        void partitionClasses(float_t threshold);

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
        void normalize(float_t border);

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


        /**
         * Destructor
         */
        virtual ~DataVector();

      private:
        /// Array to store the data
        float_t* data;

        /// Array of  correction for Kahan's summation in accumulate()
        float_t* correction;

        /// Number of elements of the data vector
        size_t size;
        /// Number of additional rows for which memory has already been reserved
        size_t unused;
        /// Number of elements by which the reserved memory is increased, if adding a row would exceed the storage reserved so far.
        size_t inc_elems;
    };

  }
}

#endif /* DATAVECTOR_HPP */
