// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/ScopedLock.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>
#include <string>

#ifndef _WIN32
#include <sys/ioctl.h>
#include <unistd.h>
#endif

namespace sgpp {
namespace base {

Printer::Printer()
    : verbose(DEFAULT_VERBOSITY),
      statusPrintingEnabled(true),
      statusLevel(0),
      indentationLevel(0),
      cursorInClearLine(true),
      lastMsgLength(0),
      lastDuration(0.0),
      lineLengthLimit(0),
      indentation(INDENTATION_LENGTH, INDENTATION_CHAR),
      stream(&std::cout) {
#ifndef _WIN32
  struct winsize w;

  if ((ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) != -1) && (w.ws_col > 0)) {
    lineLengthLimit = w.ws_col - 1;
  }

#endif
}

Printer& Printer::getInstance() {
  static Printer printer;
  return printer;
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
  (*stream) << msg;
  printStatusNewLine();

  // timing
  watches.push(base::SGppStopwatch());

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

  std::string printMsg = msg;
  const size_t totalIndentationLength = indentationLevel * INDENTATION_LENGTH;

  if ((lineLengthLimit > 0) && (totalIndentationLength + printMsg.length() > lineLengthLimit)) {
    printMsg = printMsg.substr(0, lineLengthLimit - totalIndentationLength - 3) + "...";
  }

  // print indentation
  if (cursorInClearLine) {
    printStatusIdentation();
  }

  // go back to start of last message and overwrite it with the new message
  (*stream) << std::string(lastMsgLength, '\b') << printMsg;

  // new message is too short ==> print spaces to hide the old one and
  // go back again
  if (lastMsgLength > printMsg.length()) {
    (*stream) << std::string(lastMsgLength - printMsg.length(), ' ')
              << std::string(lastMsgLength - printMsg.length(), '\b');
  }

  (*stream) << std::flush;
  lastMsgLength = printMsg.length();
  cursorInClearLine = false;
}

void Printer::printStatusNewLine() {
  if (!statusPrintingEnabled || (statusLevel > verbose)) {
    // status printing disabled or verbose level too low
    return;
  }

  (*stream) << std::endl;
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
    (*stream) << indentation;
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
  base::SGppStopwatch watch = watches.top();
  lastDuration = watch.stop();

  // print message
  std::string time_msg =
      "Done in " + std::to_string(static_cast<size_t>(1000.0 * lastDuration)) + "ms";

  if (msg == "") {
    (*stream) << time_msg << ".";
  } else {
    (*stream) << time_msg << ", " << msg;
  }

  printStatusNewLine();
  lastMsgLength = 0;
  watches.pop();
  cursorInClearLine = true;
}

void Printer::enableStatusPrinting() {
  ScopedLock lock(mutex);
  statusPrintingEnabled = true;
}

void Printer::disableStatusPrinting() {
  ScopedLock lock(mutex);
  statusPrintingEnabled = false;
}

bool Printer::isStatusPrintingEnabled() {
  ScopedLock lock(mutex);
  return statusPrintingEnabled;
}

double Printer::getLastDurationSecs() const { return lastDuration; }

MutexType& Printer::getMutex() { return mutex; }

std::ostream* Printer::getStream() const { return stream; }

void Printer::setStream(std::ostream* stream) { this->stream = stream; }

size_t Printer::getLineLengthLimit() { return lineLengthLimit; }

void Printer::setLineLengthLimit(size_t lineLengthLimit) {
  this->lineLengthLimit = lineLengthLimit;
}

void Printer::printSLE(SLE& system) const {
  const size_t n = system.getDimension();

  (*stream) << "A = [";

  // print matrix
  for (size_t i = 0; i < n; i++) {
    if (i > 0) {
      (*stream) << ";\n";
    }

    for (size_t j = 0; j < n; j++) {
      if (j > 0) {
        (*stream) << ", ";
      }

      (*stream) << system.getMatrixEntry(i, j);
    }
  }

  (*stream) << "]\n";
}
}  // namespace base
}  // namespace sgpp
