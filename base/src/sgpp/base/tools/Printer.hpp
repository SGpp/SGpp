// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/MutexType.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cstddef>
#include <stack>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * Concatenate output stream with std::vector.
 *
 * @param stream    output stream
 * @param x         vector
 * @return          stream
 */
template <class T>
inline std::ostream& operator<<(std::ostream& stream, const std::vector<T>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    stream << ((i > 0) ? ", " : "[") << x[i];
  }

  return stream << "]";
}

/**
 * Concatenate output stream with DataVector.
 *
 * @param stream    output stream
 * @param x         vector
 * @return          stream
 */
inline std::ostream& operator<<(std::ostream& stream, const DataVector& x) {
  return stream << x.toString();
}

/**
 * Concatenate output stream with GridPoint.
 *
 * @param stream    output stream
 * @param x         pointer to grid point
 * @return          stream
 */
inline std::ostream& operator<<(std::ostream& stream, const GridPoint& x) {
  DataVector xCoord(x.getDimension());

  for (size_t t = 0; t < x.getDimension(); t++) {
    xCoord[t] = x.getStandardCoordinate(t);
  }

  return stream << xCoord;
}

/**
 * Singleton class to facilitate debugging output.
 * Use with the sgpp::printer instance.
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
   * @return singleton instance
   */
  static Printer& getInstance();

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
   * @return          whether status printing is enabled or not
   */
  bool isStatusPrintingEnabled();

  /**
   * @return          current verbosity level
   */
  inline int getVerbosity() const { return verbose; }

  /**
   * @param level     new verbosity level
   */
  inline void setVerbosity(int level) { verbose = level; }

  /**
   * @return  running-time of the last operation terminated by
   *          printStatusEnd()
   */
  double getLastDurationSecs() const;

  /**
   * @return  internal mutex
   */
  MutexType& getMutex();

  /**
   * @return stream used for printing (default std::cout)
   */
  std::ostream* getStream() const;

  /**
   * @param stream stream used for printing (default std::cout)
   */
  void setStream(std::ostream* stream);

  /**
   * @return maximum length of lines, 0 if unbounded
   */
  size_t getLineLengthLimit();

  /**
   * @param lineLengthLimit maximum length of lines, 0 if unbounded
   */
  void setLineLengthLimit(size_t lineLengthLimit);

  /**
   * Print a system of linear equations.
   *
   * @param system        system to be printed
   */
  void printSLE(SLE& system) const;

 protected:
  static const size_t INDENTATION_LENGTH = 4;
  static const char INDENTATION_CHAR = ' ';

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
  /// stack of stop watches (started at time of printStatusBegin() calls)
  std::stack<SGppStopwatch> watches;
  /// length of last operation in seconds
  double lastDuration;
  /// maximum length of lines, 0 if unbounded
  size_t lineLengthLimit;
  /// indentation
  std::string indentation;

  /// internal mutex
  MutexType mutex;

  /// stream used for printing (default std::cout)
  std::ostream* stream;

 private:
  /**
   * Constructor.
   */
  Printer();

  /**
   * Deleted copy constructor.
   */
  Printer(const Printer&) = delete;

  /**
   * Deleted assignment operator.
   */
  void operator=(const Printer&) = delete;
};
}  // namespace base
}  // namespace sgpp
