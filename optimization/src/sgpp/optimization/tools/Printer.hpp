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
#include <sgpp/optimization/sle/system/SLE.hpp>
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
    inline std::ostream& operator<<(std::ostream& stream,
                                    const std::vector<T>& x) {
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
    inline std::ostream& operator<<(std::ostream& stream,
                                    const base::DataVector& x) {
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
    inline std::ostream& operator<<(std::ostream& stream,
                                    SGPP::base::GridIndex& x) {
      for (size_t t = 0; t < x.dim(); t++) {
        stream << ((t > 0) ? ", " : "[") << x.getCoord(t);
      }

      return stream << "]";
    }

    /**
     * Singleton class to facilitate debugging output.
     * Use with the SGPP::printer instance.
     *
     * The status printing functions won't print anything if
     * - status printing is disabled with disableStatusPrinting() OR
     * - the current "indentation level" (calls of printStatusBegin())
     * exceeds the "verbosity".
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
         * The last status is erased,
         * if printStatusNewLine has not been called before.
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
         * Retrieve the running time of the operation with
         * getLastDurationSecs()
         * (also works with nested calls of printStatusBegin()).
         *
         * @param msg   short description of the result,
         *              e.g. "success" or "error" (optional)
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
         * @return  running-time of the last operation terminated by
         *          printStatusEnd()
         */
        float_t getLastDurationSecs() const;

        /**
         * @return  internal mutex
         */
        MutexType& getMutex();

        /**
         * Print a grid (grid points and function values).
         *
         * @param gridGen       grid to be printed
         */
        void printIterativeGridGenerator(
          const IterativeGridGenerator& gridGen) const;

        /**
         * Print a system of linear equations.
         *
         * @param system        system to be printed
         */
        void printSLE(SLE& system) const;

      protected:
        /// verbosity level
        int verbose;
        /// whether status printing is enabled
        bool statusPrintingEnabled;

        /// current status level
        int statusLevel;
        /// current indentation level (can be less than statusLevel)
        int indentationLevel;
        /// is current line empty? (false if there's an old status update)
        bool cursorInClearLine;
        /// length of the last status message in characters
        size_t lastMsgLength;
        /// stack of the starting times (time of printStatusBegin() calls)
        std::stack<timespec> startTimes;
        /// length of last operation in seconds
        float_t lastDuration;

        /// internal mutex
        MutexType mutex;
    };

    /// singleton printer instance
    extern Printer printer;

  }
}

#endif /* SGPP_OPTIMIZATION_TOOLS_PRINTER_HPP */
