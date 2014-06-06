/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/application/ScreenOutput.hpp"

namespace sg {
  namespace base {

#ifdef _WIN32
    ScreenOutput::ScreenOutput() {
      hCon = GetStdHandle(STD_OUTPUT_HANDLE);
      first_run = true;
    }

    void ScreenOutput::update(size_t progress, std::string status) {
      int i;
      GetConsoleScreenBufferInfo(hCon, &info);

      if (!first_run) {
        pos.Y = info.dwCursorPosition.Y - 3;
        pos.X = 0;
        SetConsoleCursorPosition(hCon, pos);
      } else {
        first_run = false;
      }

      std::cout << "[";

      for (i = 0; i < ((int)(((double)progress) * 0.64)); i++) {
        std::cout << '\xdb';
      }

      for (; i < 64; i++) {
        std::cout << " ";
      }

      std::cout << "]  " << progress << "%" << std::endl << std::endl;
      std::cout << status << "       " << std::endl;
    }
#endif

#ifndef _WIN32
    ScreenOutput::ScreenOutput() {
      first_run = true;
    }

    void ScreenOutput::update(size_t progress, std::string status) {
      size_t i;

      if (!first_run) {
        std::cout << "\033[4A" << std::endl;
      } else {
        first_run = false;
      }

      std::cout << "[";

      for (i = 0; i < ((int)(((double)progress) * 0.64)); i++) {
        std::cout << "=";
      }

      for (; i < 64; i++) {
        std::cout << " ";
      }

      std::cout << "]  " << progress << "%" << std::endl << std::endl;
      std::cout << status << "          " << std::endl;
    }
#endif

    ScreenOutput::~ScreenOutput() {
    }

    void ScreenOutput::writeTitle(std::string appTitle, std::string appAuthor) {
      std::cout << std::endl << std::endl;

      std::string empty = "                                                      ";
      appTitle.append(empty);
      appAuthor.append(empty);

      appTitle = appTitle.substr(0, empty.length() - 2);
      appAuthor = appAuthor.substr(0, empty.length());

      std::cout << appTitle << "########  ##########" << std::endl;
      std::cout << "=========================================           " << "  ##  ##  ##  ##  ##" << std::endl;
      std::cout << "                                                      ##  ##  ##  ##  ##" << std::endl;
      std::cout << appAuthor << "##  ##  ##  ##  ##" << std::endl;
      std::cout << "                                                      ##  ######  ##  ##" << std::endl << std::endl;
    }

    void ScreenOutput::writeHelp(std::string helpText) {
      std::cout << helpText << std::endl << std::endl;
    }

    void ScreenOutput::writeStartSolve(std::string text) {
      std::string line;
      size_t line_length = (size_t)text.length();

      while (line_length > 0) {
        line.append("-");
        line_length--;
      }

      std::cout << std::endl;
      std::cout << text << std::endl;
      std::cout << line << std::endl;
      std::cout << std::endl;
    }

    void ScreenOutput::writeEmptyLines(size_t numLines) {
      for (size_t i = 0; i < numLines; i++) {
        std::cout << std::endl;
      }
    }

  }
}
