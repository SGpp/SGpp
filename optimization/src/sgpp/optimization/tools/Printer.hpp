// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TOOLS_PRINTER_HPP
#define SGPP_OPTIMIZATION_TOOLS_PRINTER_HPP

#include <algorithm>
#include <cstddef>
#include <stack>
#include <sys/time.h>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>
#include <sgpp/optimization/sle/system/System.hpp>
#include <sgpp/optimization/tools/MutexType.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Concatenate output stream with std::vector.
     *
     * @param stream    output stream
     * @param x         vector
     */
    template <class T>
    inline std::ostream& operator<<(std::ostream& stream, const std::vector<T>& x) {
      for (size_t i = 0; i < x.size(); i++) {
        stream << ((i > 0) ? ", " : "[") << x[i];
      }

      return stream << "]";
    }

    /**
     * Concatenate output stream with base::DataVector.
     *
     * @param stream    output stream
     * @param x         vector
     */
    inline std::ostream& operator<<(std::ostream& stream, const base::DataVector& x) {
      for (size_t i = 0; i < x.getSize(); i++) {
        stream << ((i > 0) ? ", " : "[") << x.get(i);
      }

      return stream << "]";
    }

    /**
     * Concatenate output stream with base::GridIndex.
     *
     * @param stream    output stream
     * @param x         pointer to grid point
     */
    inline std::ostream& operator<<(std::ostream& stream, SGPP::base::GridIndex* x) {
      for (size_t t = 0; t < x->dim(); t++) {
        stream << ((t > 0) ? ", " : "[") << x->getCoord(t);
      }

      return stream << "]";
    }

    /**
     * Convert numbers to strings.
     * (In C++11, you would use std::to_string.)
     *
     * @param x     number
     * @return      string representing x
     */
    template <typename T>
    std::string toString(const T& x) {
      std::ostringstream oss;
      oss << x;
      return oss.str();
    }

    namespace tools {

      /**
       * Singleton class to facilitate debugging output.
       * Use with the SGPP::tools::printer instance.
       *
       * The status printing functions won't print anything if
       * - status printing is disabled with disableStatusPrinting() OR
       * - the current "indentation level" (calls of printStatusBegin()) exceeds the "verbosity".
       */
      class Printer {
        public:
          /// default verbosity
          static const int DEFAULT_VERBOSITY = 0;

          /**
           * Constructor.
           */
          Printer();

          /**
           * Call at the beginning of a time-consuming operation.
           * Locks and unlocks an OpenMP mutex.
           *
           * @param msg   short description of the operation
           */
          void printStatusBegin(const std::string& msg);

          /**
           * Call for printing status updates on the operation.
           * The last status is erased, if printStatusNewLine has not been called before.
           * Locks and unlocks an OpenMP mutex.
           *
           * @param msg   status message
           */
          void printStatusUpdate(const std::string& msg);

          /**
           * End the last status update and place the cursor in a newline.
           */
          void printStatusNewLine();

          /**
           * Internal function printing the indentation.
           */
          void printStatusIdentation();

          /**
           * Call at the end of a time-consuming operation.
           * Locks and unlocks an OpenMP mutex.
           *
           * Retrieve the running time of the operation with getLastDurationSecs()
           * (also works with nested calls of printStatusBegin()).
           *
           * @param msg   short description of the result, e.g. "success" or "error" (optional)
           */
          void printStatusEnd(const std::string& msg = "");

          /**
           * Enable the printStatus... functions.
           * (They're enabled by default.)
           */
          void enableStatusPrinting();

          /**
           * Disable the printStatus... functions.
           * (They're enabled by default.)
           */
          void disableStatusPrinting();

          /**
           * Output a grid (grid points and function values) to a file.
           *
           * The format is as follows:
           *
           * <pre>
           * size_t   N (number of grid points)
           * size_t   d (dimension)
           * for j = 0, ..., N-1
           *     for t = 0, ..., d-1
           *         float_t           grid_point[j].getCoord(t)
           *         unsigned int     grid_point[j].level(t)
           *         unsigned int     grid_point[j].index(t)
           *     end
           *     float_t   function_value[j]
           * end
           * </pre>
           *
           * @param filename      filename of the file to be written
           * @param gridGen       iterative grid generator containing the grid points and function values
           */
          void printGridToFile(const std::string& filename,
                               const gridgen::IterativeGridGenerator& gridGen) const;

          /**
           * Output a base::DataVector to a file.
           * It's printMatrixToFile with the vector as one column.
           *
           * @param filename      filename of the file to be written
           * @param x             vector
           */
          void printVectorToFile(const std::string& filename, base::DataVector& x) const;

          /**
           * Output a std::vector to a file.
           * It's printMatrixToFile with the vector as one column.
           *
           * @param filename      filename of the file to be written
           * @param x             vector
           */
          template <class T>
          void printVectorToFile(const std::string& filename, const std::vector<T>& x) const {
            printMatrixToFile(filename, x, x.size(), 1);
          }

          /**
           * Output a base::DataMatrix to a file.
           *
           * @param filename      filename of the file to be written
           * @param A             matrix
           */
          void printMatrixToFile(const std::string& filename, base::DataMatrix& A) const;

          /**
           * Output a matrix (stored row-wise in a std::vector) to a file.
           *
           * The format is as follows:
           *
           * <pre>
           * size_t       m
           * size_t       n
           * char[16]     type string (one of "uint8", "uint16", "uint32", "uint64", "float_t", "string"
           *              or "other", right-padded with spaces to 16 characters)
           * for i = 0, ..., m*n - 1
           *     T        A[i] (size depending on template parameter,
           *              strings are written null-terminatedly)
           * end
           * </pre>
           *
           * @param filename      filename of the file to be written
           * @param A             matrix
           * @param m             number of rows
           * @param n             number of columns
           */
          template <class T>
          void printMatrixToFile(const std::string& filename, const std::vector<T>& A,
                                 size_t m, size_t n) const {
            std::ofstream f(filename.c_str(), std::ios::binary);
            const char* type = getTypeString(A);

            // header (size and type)
            f.write(reinterpret_cast<const char*>(&m), sizeof(m));
            f.write(reinterpret_cast<const char*>(&n), sizeof(n));
            f.write(type, 16);

            // entries
            for (size_t i = 0; i < m * n; i++) {
              writeEntryToFile(f, A[i]);
            }

            f.close();
          }

          /**
           * Output a matrix (stored as a vector of row vectors) to a file.
           * Copies the matrix and calls the other version with the vectorized matrix.
           *
           * @param filename      filename of the file to be written
           * @param A             matrix
           */
          template <class T>
          void printMatrixToFile(const std::string& filename,
                                 const std::vector<std::vector<T> >& A) const {
            size_t m = A.size();
            size_t n = (A.empty() ? 0 : A[0].size());
            std::vector<T> B(m * n);

            for (size_t i = 0; i < m; i++) {
              std::copy(A[i].begin(), A[i].end(), B.begin() + i * n);
            }

            printMatrixToFile(filename, B, m, n);
          }

          /**
           * Print a grid (grid points and function values).
           *
           * @param gridGen       grid to be printed
           */
          void printIterativeGridGenerator(const gridgen::IterativeGridGenerator& gridGen) const;

          /**
           * Print a system of linear equation.
           *
           * @param system        system to be printed
           */
          void printSLE(sle::system::System& system) const;

          /**
           * @return          current verbosity level
           */
          inline int getVerbosity() const {
            return verbose;
          }

          /**
           * @param level     new verbosity level
           */
          inline void setVerbosity(int level) {
            verbose = level;
          }

          /**
           * @return  running-time of the last operation terminated by printStatusEnd()
           */
          float_t getLastDurationSecs() const;

          /**
           * @return  internal mutex
           */
          MutexType& getMutex();

        protected:
          /// verbosity level
          int verbose;
          /// whether status printing is enabled
          bool statusPrintingEnabled;

          /// current status level
          int statusLevel;
          /// current indentation level (can be less than status_level)
          int indentationLevel;
          /// whether the current line is clear (false if there is an old status update)
          bool cursorInClearLine;
          /// length of the last status message in characters
          size_t lastMsgLength;
          /// stack of the starting times of the operations (time points of printStatusBegin() calls)
          std::stack<timespec> startTimes;
          /// length of last operation in seconds
          float_t lastDuration;

          /// internal mutex
          MutexType mutex;

          /**
           * @return  type string for unknown types (right-padded "other")
           */
          template <class T>
          const char* getTypeString(const std::vector<T>& A) const {
            (void)A;
            return "other           ";
          }

          /**
           * Write bytes of entry representation to output stream (called by printMatrixToFile()).
           *
           * @param f         output stream to be written to
           * @param entry     entry to be written
           */
          template <class T>
          void writeEntryToFile(std::ofstream& f, const T& entry) const {
            f.write(reinterpret_cast<const char*>(&entry), sizeof(T));
          }
      };

      /**
       * @return  type string for uint8_t (right-padded "uint8")
       */
      template <>
      const char* Printer::getTypeString(const std::vector<uint8_t>& A) const;

      /**
       * @return  type string for uint16_t (right-padded "uint16")
       */
      template <>
      const char* Printer::getTypeString(const std::vector<uint16_t>& A) const;

      /**
       * @return  type string for uint32_t (right-padded "uint32")
       */
      template <>
      const char* Printer::getTypeString(const std::vector<uint32_t>& A) const;

      /**
       * @return  type string for uint64_t (right-padded "uint64")
       */
      template <>
      const char* Printer::getTypeString(const std::vector<uint64_t>& A) const;

      /**
       * @return  type string for float_t (right-padded "float_t")
       */
      template <>
      const char* Printer::getTypeString(const std::vector<float_t>& A) const;

      /**
       * @return  type string for string (right-padded "string")
       */
      template <>
      const char* Printer::getTypeString(const std::vector<std::string>& A) const;

      /**
       * Write string to output stream (called by printMatrixToFile()).
       *
       * @param f         output stream to be written to
       * @param entry     entry to be written
       */
      template <>
      void Printer::writeEntryToFile(std::ofstream& f, const std::string& entry) const;

      /// singleton printer instance
      extern Printer printer;

    }
  }
}

#endif
