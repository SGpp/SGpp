// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/ScopedLock.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <iostream>

namespace SGPP {
  namespace optimization {

    Printer printer;

    Printer::Printer() :
      verbose(DEFAULT_VERBOSITY),
      statusPrintingEnabled(true),
      statusLevel(0),
      indentationLevel(0),
      cursorInClearLine(true),
      lastMsgLength(0),
      lastDuration(0.0) {
    }

    void Printer::printStatusBegin(const std::string& msg) {
      ScopedLock lock(mutex);

      if (!statusPrintingEnabled || (statusLevel > verbose)) {
        // status printing disabled or verbose level too low
        statusLevel++;
        return;
      }

      // go to new line
      if (!cursorInClearLine) {
        printStatusNewLine();
      }

      // print indentation
      printStatusIdentation();

      // print message
      std::cout << msg;
      printStatusNewLine();

      // timing
      timespec ts;
      clock_gettime(CLOCK_MONOTONIC, &ts);
      startTimes.push(ts);

      lastMsgLength = 0;
      cursorInClearLine = true;
      statusLevel++;
      indentationLevel++;
    }

    void Printer::printStatusUpdate(const std::string& msg) {
      ScopedLock lock(mutex);

      if (!statusPrintingEnabled || (statusLevel > verbose)) {
        // status printing disabled or verbose level too low
        return;
      }

      // print indentation
      if (cursorInClearLine) {
        printStatusIdentation();
      }

      // go back to start of last message and overwrite it with the new message
      std::cout << std::string(lastMsgLength, '\b') << msg;

      // new message is too short ==> print spaces to hide the old one and
      // go back again
      if (lastMsgLength > msg.length()) {
        std::cout << std::string(lastMsgLength - msg.length(), ' ') <<
                  std::string(lastMsgLength - msg.length(), '\b');
      }

      std::cout << std::flush;
      lastMsgLength = msg.length();
      cursorInClearLine = false;
    }

    void Printer::printStatusNewLine() {
      if (!statusPrintingEnabled || (statusLevel > verbose)) {
        // status printing disabled or verbose level too low
        return;
      }

      std::cout << std::endl;
      lastMsgLength = 0;
      cursorInClearLine = true;
    }

    void Printer::printStatusIdentation() {
      if (!statusPrintingEnabled || (statusLevel > verbose)) {
        // status printing disabled or verbose level too low
        return;
      }

      // print indentation
      for (int i = 0; i < indentationLevel; i++) {
        std::cout << "    ";
      }
    }

    void Printer::printStatusEnd(const std::string& msg) {
      ScopedLock lock(mutex);

      statusLevel--;

      if (!statusPrintingEnabled || (statusLevel > verbose)) {
        // status printing disabled or verbose level too low
        return;
      }

      // go to new line
      if (!cursorInClearLine) {
        printStatusNewLine();
      }

      // print indentation
      indentationLevel--;
      printStatusIdentation();

      // timing
      timespec ts1 = startTimes.top(), ts2;
      clock_gettime(CLOCK_MONOTONIC, &ts2);
      lastDuration = static_cast<float_t>(ts2.tv_sec - ts1.tv_sec) +
                     static_cast<float_t>(ts2.tv_nsec - ts1.tv_nsec) / 1e9;

      // print message
      std::string time_msg =
        "Done in " +
        std::to_string(static_cast<size_t>(1000.0 * lastDuration)) + "ms";

      if (msg == "") {
        std::cout << time_msg << ".";
      } else {
        std::cout << time_msg << ", " << msg;
      }

      printStatusNewLine();
      lastMsgLength = 0;
      startTimes.pop();
      cursorInClearLine = true;
    }

    void Printer::enableStatusPrinting() {
      statusPrintingEnabled = true;
    }

    void Printer::disableStatusPrinting() {
      statusPrintingEnabled = false;
    }

    bool Printer::isStatusPrintingEnabled() {
      return statusPrintingEnabled;
    }

    float_t Printer::getLastDurationSecs() const {
      return lastDuration;
    }

    MutexType& Printer::getMutex() {
      return mutex;
    }

    void Printer::printIterativeGridGenerator(
      const IterativeGridGenerator& grid_gen) const {
      base::GridStorage& gridStorage = *grid_gen.getGrid().getStorage();
      const base::DataVector& functionValues =
        grid_gen.getFunctionValues();

      for (size_t i = 0; i < gridStorage.size(); i++) {
        if (i > 0) {
          std::cout << "\n";
        }

        // print grid point and function value
        std::cout << *gridStorage.get(i) << ", " << functionValues.get(i);
      }
    }

    void Printer::printSLE(SLE& system) const {
      const size_t n = system.getDimension();

      std::cout << "A = [";

      // print matrix
      for (size_t i = 0; i < n; i++) {
        if (i > 0) {
          std::cout << ";\n";
        }

        for (size_t j = 0; j < n; j++) {
          if (j > 0) {
            std::cout << ", ";
          }

          std::cout << system.getMatrixEntry(i, j);
        }
      }

      std::cout << "]";
    }

  }
}
