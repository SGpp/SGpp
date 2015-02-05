// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/optimization/tools/ScopedLock.hpp>

#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace tools {

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

        // new message is too short ==> print spaces to hide the old one and go back again
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
          "Done in " + toString(static_cast<size_t>(1000.0 * lastDuration)) + "ms";

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

      void Printer::printGridToFile(const std::string& filename,
                                    const gridgen::IterativeGridGenerator& gridGen) const {
        base::GridStorage* grid_storage = gridGen.getGrid().getStorage();
        const std::vector<float_t>& function_values = gridGen.getFunctionValues();
        size_t N = grid_storage->size();
        size_t d = grid_storage->dim();

        std::ofstream f(filename.c_str(), std::ios::out | std::ios::binary);

        // header (number of points and dimensions)
        f.write(reinterpret_cast<const char*>(&N), sizeof(N));
        f.write(reinterpret_cast<const char*>(&d), sizeof(d));

        for (size_t j = 0; j < N; j++) {
          base::GridIndex* gp = grid_storage->get(j);

          for (size_t t = 0; t < d; t++) {
            float_t x = gp->getCoord(t);
            base::GridIndex::level_type l = gp->getLevel(t);
            base::GridIndex::index_type i = gp->getIndex(t);

            // coordinate, level and index of current grid point
            f.write(reinterpret_cast<const char*>(&x), sizeof(x));
            f.write(reinterpret_cast<const char*>(&l), sizeof(l));
            f.write(reinterpret_cast<const char*>(&i), sizeof(i));
          }

          // function value at the current grid point
          f.write(reinterpret_cast<const char*>(&function_values[j]), sizeof(float_t));
        }

        f.close();
      }

      void Printer::printVectorToFile(const std::string& filename, base::DataVector& x) const {
        // convert DataVector to std::vector
        std::vector<float_t> xVector(x.getPointer(), x.getPointer() + x.getSize());
        printMatrixToFile(filename, xVector, x.getSize(), 1);
      }

      void Printer::printMatrixToFile(const std::string& filename, base::DataMatrix& A) const {
        // convert DataMatrix to std::vector
        std::vector<float_t> AVector(A.getPointer(), A.getPointer() + A.getNrows() * A.getNcols());
        printMatrixToFile(filename, AVector, A.getNrows(), A.getNcols());
      }

      void Printer::printIterativeGridGenerator(const gridgen::IterativeGridGenerator& grid_gen) const {
        base::GridStorage* gridStorage = grid_gen.getGrid().getStorage();
        const std::vector<float_t>& functionValues = grid_gen.getFunctionValues();

        for (size_t i = 0; i < gridStorage->size(); i++) {
          if (i > 0) {
            std::cout << "\n";
          }

          // print grid point and function value
          std::cout << gridStorage->get(i) << ", " << functionValues[i];
        }
      }

      void Printer::printSLE(sle::system::System& system) const {
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

      float_t Printer::getLastDurationSecs() const {
        return lastDuration;
      }

      MutexType& Printer::getMutex() {
        return mutex;
      }

      template <>
      const char* Printer::getTypeString(const std::vector<uint8_t>& A) const {
        (void)A;
        return "uint8           ";
      }

      template <>
      const char* Printer::getTypeString(const std::vector<uint16_t>& A) const {
        (void)A;
        return "uint16          ";
      }

      template <>
      const char* Printer::getTypeString(const std::vector<uint32_t>& A) const {
        (void)A;
        return "uint32          ";
      }

      template <>
      const char* Printer::getTypeString(const std::vector<uint64_t>& A) const {
        (void)A;
        return "uint64          ";
      }

      template <>
      const char* Printer::getTypeString(const std::vector<float_t>& A) const {
        (void)A;
        return "float_t          ";
      }

      template <>
      const char* Printer::getTypeString(const std::vector<std::string>& A) const {
        (void)A;
        return "string          ";
      }

      template <>
      void Printer::writeEntryToFile(std::ofstream& f, const std::string& entry) const {
        // write string terminated by null character
        const char null_char[1] = {'\0'};
        f << entry;
        f.write(null_char, 1);
      }

    }
  }
}
